//! dbNSFP parser for building .osa annotation files.
//!
//! dbNSFP provides pre-computed functional predictions (SIFT, PolyPhen,
//! AlphaMissense, ESM1b, REVEL, CADD, etc.) for all possible missense variants.
//!
//! This parser extracts:
//!   * SIFT, PolyPhen2 HDIV  -- existing fields (single-value, "first ;-split")
//!   * AlphaMissense + class -- per-transcript; collapsed by MAX across isoforms
//!   * ESM1b + class         -- per-transcript LLR; collapsed by MIN (most damaging)
//!   * REVEL                 -- per-variant score (single value)
//!
//! Rationale for the collapse rule: AlphaMissense and ESM1b are isoform-specific
//! (Brandes 2023 reports ~2M variants damaging only in some isoforms). Taking the
//! worst score across transcripts is the conservative pathogenicity reading, and
//! matches the `max_alphamissense_any_tx` / `min_esm1b_any_tx` secondary columns
//! used downstream by R/tier_variants.R. SIFT/PolyPhen/REVEL retain "first" for
//! backward compatibility.

use crate::common::AnnotationRecord;
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::BufRead;

/// Parse a dbNSFP TSV file to extract SIFT, PolyPhen, AlphaMissense, ESM1b and REVEL.
///
/// Expected header columns (any subset; columns absent in this dbNSFP build are skipped):
/// `#chr`, `pos(1-based)`, `ref`, `alt`,
/// `SIFT_score`, `SIFT_pred`,
/// `Polyphen2_HDIV_score`, `Polyphen2_HDIV_pred`,
/// `AlphaMissense_score`, `AlphaMissense_pred`,
/// `ESM1b_score`, `ESM1b_pred`,
/// `REVEL_score`
///
/// Column indices are auto-detected from the header row (case-insensitive).
pub fn parse_dbnsfp<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
) -> Result<Vec<AnnotationRecord>> {
    let mut records = Vec::new();
    let mut col_indices: Option<DbNsfpColumns> = None;

    for line in reader.lines() {
        let line = line.context("Reading dbNSFP line")?;

        if line.starts_with('#') || line.starts_with("chr\t") {
            // Parse header to find column indices
            let header = line.trim_start_matches('#');
            col_indices = Some(DbNsfpColumns::from_header(header)?);
            continue;
        }

        if line.is_empty() {
            continue;
        }

        let cols = match &col_indices {
            Some(c) => c,
            None => continue,
        };

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() <= cols.max_idx() {
            continue;
        }

        let chrom = normalize_chrom(fields[cols.chr]);
        let chrom_idx = match chrom_to_idx.get(&chrom) {
            Some(&idx) => idx,
            None => continue,
        };

        let pos: u32 = match fields[cols.pos].parse() {
            Ok(p) => p,
            Err(_) => continue,
        };

        let ref_allele = fields[cols.ref_col].to_string();
        let alt_allele = fields[cols.alt].to_string();

        let mut parts = Vec::new();

        // SIFT (first ;-split value, backward-compatible)
        if let Some(idx) = cols.sift_score {
            if let Some((score, pred)) = first_score_pred(&fields, idx, cols.sift_pred) {
                let pred_str = match pred.as_deref() {
                    Some("D") => "deleterious",
                    Some("T") => "tolerated",
                    _ => "",
                };
                if !pred_str.is_empty() {
                    parts.push(format!("\"sift\":\"{}({:.3})\"", pred_str, score));
                } else {
                    parts.push(format!("\"sift\":\"{:.3}\"", score));
                }
            }
        }

        // PolyPhen2 HDIV (first ;-split value, backward-compatible)
        if let Some(idx) = cols.polyphen_score {
            if let Some((score, pred)) = first_score_pred(&fields, idx, cols.polyphen_pred) {
                let pred_str = match pred.as_deref() {
                    Some("D") => "probably_damaging",
                    Some("P") => "possibly_damaging",
                    Some("B") => "benign",
                    _ => "",
                };
                if !pred_str.is_empty() {
                    parts.push(format!("\"polyphen\":\"{}({:.3})\"", pred_str, score));
                } else {
                    parts.push(format!("\"polyphen\":\"{:.3}\"", score));
                }
            }
        }

        // AlphaMissense -- MAX across per-transcript ;-split values (worst case).
        // dbNSFP AlphaMissense_pred letters: P=likely_pathogenic, A=ambiguous, B=likely_benign.
        if let Some(idx) = cols.alphamissense_score {
            if let Some((score, pred)) =
                worst_score_pred(&fields, idx, cols.alphamissense_pred, WorstDir::Max)
            {
                let class = match pred.as_deref() {
                    Some("P") => "likely_pathogenic",
                    Some("A") => "ambiguous",
                    Some("B") => "likely_benign",
                    _ => "",
                };
                parts.push(format!("\"alphamissense\":{:.4}", score));
                if !class.is_empty() {
                    parts.push(format!("\"am_class\":\"{}\"", class));
                }
            }
        }

        // ESM1b -- MIN across per-transcript ;-split values (most damaging LLR).
        // dbNSFP ESM1b_pred letters: D=damaging, T=tolerated.
        if let Some(idx) = cols.esm1b_score {
            if let Some((score, pred)) =
                worst_score_pred(&fields, idx, cols.esm1b_pred, WorstDir::Min)
            {
                let class = match pred.as_deref() {
                    Some("D") => "damaging",
                    Some("T") => "tolerated",
                    _ => "",
                };
                parts.push(format!("\"esm1b\":{:.4}", score));
                if !class.is_empty() {
                    parts.push(format!("\"esm1b_class\":\"{}\"", class));
                }
            }
        }

        // REVEL -- single per-variant score; first ;-split value if multi-cardinality.
        if let Some(idx) = cols.revel_score {
            if let Some((score, _)) = first_score_pred(&fields, idx, None) {
                parts.push(format!("\"revel\":{:.4}", score));
            }
        }

        if parts.is_empty() {
            continue;
        }

        records.push(AnnotationRecord {
            chrom_idx,
            position: pos,
            ref_allele,
            alt_allele,
            json: format!("{{{}}}", parts.join(",")),
        });
    }

    records.sort_by(|a, b| a.chrom_idx.cmp(&b.chrom_idx).then(a.position.cmp(&b.position)));
    Ok(records)
}

/// Per-field collapse direction for per-transcript `;`-split cells.
#[derive(Clone, Copy, PartialEq, Eq)]
enum WorstDir {
    Max,
    Min,
}

/// Read the first non-`.` `;`-split element of a score cell and its paired pred letter.
fn first_score_pred(
    fields: &[&str],
    score_idx: usize,
    pred_idx: Option<usize>,
) -> Option<(f64, Option<String>)> {
    let score_cell = fields.get(score_idx)?;
    if *score_cell == "." {
        return None;
    }
    let score: f64 = score_cell.split(';').next()?.parse().ok()?;
    let pred = pred_idx.and_then(|i| {
        let cell = fields.get(i)?;
        if *cell == "." {
            return None;
        }
        cell.split(';').next().and_then(|p| {
            if p == "." || p.is_empty() {
                None
            } else {
                Some(p.to_string())
            }
        })
    });
    Some((score, pred))
}

/// Read all `;`-split elements of a score cell and return the MIN or MAX score,
/// together with the paired pred letter at that same index (best-effort
/// alignment with the pred column, which may have a different cardinality).
fn worst_score_pred(
    fields: &[&str],
    score_idx: usize,
    pred_idx: Option<usize>,
    dir: WorstDir,
) -> Option<(f64, Option<String>)> {
    let score_cell = fields.get(score_idx)?;
    if *score_cell == "." {
        return None;
    }
    let preds: Vec<&str> = pred_idx
        .and_then(|i| fields.get(i).copied())
        .filter(|c| *c != ".")
        .map(|c| c.split(';').collect::<Vec<_>>())
        .unwrap_or_default();

    let mut best: Option<(usize, f64)> = None;
    for (i, raw) in score_cell.split(';').enumerate() {
        let raw = raw.trim();
        if raw == "." || raw.is_empty() {
            continue;
        }
        let v: f64 = match raw.parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let take = match (dir, &best) {
            (_, None) => true,
            (WorstDir::Max, Some((_, b))) => v > *b,
            (WorstDir::Min, Some((_, b))) => v < *b,
        };
        if take {
            best = Some((i, v));
        }
    }
    let (best_idx, best_score) = best?;
    let pred = preds.get(best_idx).and_then(|p| {
        let p = p.trim();
        if p == "." || p.is_empty() {
            None
        } else {
            Some(p.to_string())
        }
    });
    Some((best_score, pred))
}

struct DbNsfpColumns {
    chr: usize,
    pos: usize,
    ref_col: usize,
    alt: usize,
    sift_score: Option<usize>,
    sift_pred: Option<usize>,
    polyphen_score: Option<usize>,
    polyphen_pred: Option<usize>,
    alphamissense_score: Option<usize>,
    alphamissense_pred: Option<usize>,
    esm1b_score: Option<usize>,
    esm1b_pred: Option<usize>,
    revel_score: Option<usize>,
}

impl DbNsfpColumns {
    fn from_header(header: &str) -> Result<Self> {
        let fields: Vec<&str> = header.split('\t').collect();
        let find = |names: &[&str]| -> Option<usize> {
            fields.iter().position(|f| {
                let fl = f.to_lowercase();
                names.iter().any(|n| fl == *n)
            })
        };

        Ok(Self {
            chr: find(&["chr", "#chr"]).unwrap_or(0),
            pos: find(&["pos(1-based)", "pos", "hg38_pos"]).unwrap_or(1),
            ref_col: find(&["ref", "ref_allele"]).unwrap_or(2),
            alt: find(&["alt", "alt_allele"]).unwrap_or(3),
            sift_score: find(&["sift_score"]),
            sift_pred: find(&["sift_pred"]),
            polyphen_score: find(&["polyphen2_hdiv_score"]),
            polyphen_pred: find(&["polyphen2_hdiv_pred"]),
            alphamissense_score: find(&["alphamissense_score"]),
            alphamissense_pred: find(&["alphamissense_pred"]),
            esm1b_score: find(&["esm1b_score"]),
            esm1b_pred: find(&["esm1b_pred"]),
            revel_score: find(&["revel_score"]),
        })
    }

    fn max_idx(&self) -> usize {
        let mut m = self.alt;
        for opt in [
            self.sift_score,
            self.sift_pred,
            self.polyphen_score,
            self.polyphen_pred,
            self.alphamissense_score,
            self.alphamissense_pred,
            self.esm1b_score,
            self.esm1b_pred,
            self.revel_score,
        ] {
            if let Some(i) = opt {
                m = m.max(i);
            }
        }
        m
    }
}

fn normalize_chrom(chrom: &str) -> String {
    if chrom.starts_with("chr") {
        chrom.to_string()
    } else {
        format!("chr{}", chrom)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_dbnsfp() {
        let data = "\
#chr\tpos(1-based)\tref\talt\tSIFT_score\tSIFT_pred\tPolyphen2_HDIV_score\tPolyphen2_HDIV_pred
1\t10001\tA\tG\t0.032\tD\t0.998\tD
1\t10001\tA\tC\t0.450\tT\t0.100\tB
1\t10002\tC\tT\t.\t.\t.\t.
";
        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr1".into(), 0u16);

        let records = parse_dbnsfp(data.as_bytes(), &chrom_map).unwrap();
        // Third line has all dots, should be skipped
        assert_eq!(records.len(), 2);

        assert!(records[0].json.contains("deleterious(0.032)"));
        assert!(records[0].json.contains("probably_damaging(0.998)"));

        assert!(records[1].json.contains("tolerated(0.450)"));
        assert!(records[1].json.contains("benign(0.100)"));
    }

    #[test]
    fn test_parse_dbnsfp_extended_scores() {
        // AlphaMissense + ESM1b + REVEL, single per-transcript values.
        let data = "\
#chr\tpos(1-based)\tref\talt\tSIFT_score\tSIFT_pred\tPolyphen2_HDIV_score\tPolyphen2_HDIV_pred\tAlphaMissense_score\tAlphaMissense_pred\tESM1b_score\tESM1b_pred\tREVEL_score
17\t43045712\tC\tT\t0.001\tD\t0.999\tD\t0.9912\tP\t-12.4\tD\t0.842
22\t29683017\tA\tG\t0.150\tT\t0.300\tP\t0.6100\tA\t-6.1\tT\t0.520
";
        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr17".into(), 16);
        chrom_map.insert("chr22".into(), 21);

        let records = parse_dbnsfp(data.as_bytes(), &chrom_map).unwrap();
        assert_eq!(records.len(), 2);

        let first: serde_json::Value = serde_json::from_str(&records[0].json).unwrap();
        assert_eq!(first["alphamissense"], 0.9912);
        assert_eq!(first["am_class"], "likely_pathogenic");
        assert_eq!(first["esm1b"], -12.4);
        assert_eq!(first["esm1b_class"], "damaging");
        assert_eq!(first["revel"], 0.842);
        // Backward-compatible fields still present.
        assert!(records[0].json.contains("deleterious(0.001)"));
        assert!(records[0].json.contains("probably_damaging(0.999)"));

        let second: serde_json::Value = serde_json::from_str(&records[1].json).unwrap();
        assert_eq!(second["am_class"], "ambiguous");
        assert_eq!(second["esm1b_class"], "tolerated");
    }

    #[test]
    fn test_alphamissense_collapses_to_max_across_transcripts() {
        // ;-split: AlphaMissense should pick MAX (worst-case pathogenicity), ESM1b MIN.
        // Tied pred indices: 2nd transcript is the pathogenic one for AM (0.95);
        // 1st transcript is the damaging one for ESM1b (-12.4).
        let data = "\
#chr\tpos(1-based)\tref\talt\tAlphaMissense_score\tAlphaMissense_pred\tESM1b_score\tESM1b_pred
17\t43045712\tC\tT\t0.41;0.95\tA;P\t-12.4;-2.1\tD;T
";
        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr17".into(), 16);

        let records = parse_dbnsfp(data.as_bytes(), &chrom_map).unwrap();
        assert_eq!(records.len(), 1);

        let v: serde_json::Value = serde_json::from_str(&records[0].json).unwrap();
        // AM: max(0.41, 0.95) = 0.95, class from 2nd position = "P" -> likely_pathogenic
        assert_eq!(v["alphamissense"], 0.95);
        assert_eq!(v["am_class"], "likely_pathogenic");
        // ESM1b: min(-12.4, -2.1) = -12.4, class from 1st position = "D" -> damaging
        assert_eq!(v["esm1b"], -12.4);
        assert_eq!(v["esm1b_class"], "damaging");
    }

    #[test]
    fn test_dbnsfp_handles_missing_per_transcript_cell() {
        // First transcript dot, second has score. Collapse should pick the
        // non-missing one, no panic, no spurious zero.
        let data = "\
#chr\tpos(1-based)\tref\talt\tAlphaMissense_score\tAlphaMissense_pred\tESM1b_score
17\t43045712\tC\tT\t.;0.71\t.;A\t.;-8.1
";
        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr17".into(), 16);

        let records = parse_dbnsfp(data.as_bytes(), &chrom_map).unwrap();
        assert_eq!(records.len(), 1);

        let v: serde_json::Value = serde_json::from_str(&records[0].json).unwrap();
        assert_eq!(v["alphamissense"], 0.71);
        assert_eq!(v["am_class"], "ambiguous");
        assert_eq!(v["esm1b"], -8.1);
    }
}
