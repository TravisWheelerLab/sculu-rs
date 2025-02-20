use crate::RmBlastOutput;
use petgraph::{
    graph::UnGraph,
    unionfind::UnionFind,
    visit::{
        EdgeRef, IntoEdgeReferences, IntoNodeReferences, NodeCompactIndexable, NodeRef,
    },
};
use std::{
    cmp::{max, min},
    collections::HashMap,
};

// Make parameters?
const MIN_LEN_SIMILARITY: f64 = 0.9;
const MIN_ALIGN_COVER: f64 = 0.9;

// --------------------------------------------------
/// Given a slice of `RmBlastOutput` records, get the connected components
/// of a graph where edges are created between sequences that have been aligned.
///
/// The components are returned as `Vec<Vec<&str>>`, where each inner `Vec`
/// contains the names of the sequences in the component.
pub fn connected_components(records: Vec<RmBlastOutput>) -> Vec<Vec<String>> {
    let mut seq_to_id: HashMap<String, usize> = HashMap::new();
    let mut id_to_seq: HashMap<usize, String> = HashMap::new();
    let mut id_cnt = 0usize;

    records.iter().for_each(|record| {
        for s in &[&record.query, &record.target] {
            seq_to_id.entry(s.to_string()).or_insert_with(|| {
                let id = id_cnt;
                id_to_seq.insert(id, s.to_string());
                id_cnt += 1;
                id
            });
        }
    });


    let edges: Vec<_> = records
        .iter()
        .filter(|record| record.query != record.target)
        .filter_map(|record| {
            // The coverage of the hits to either the query or subject must
            // be at least half of the sequence length.
            let query_span = 1 + record.query_end.abs_diff(record.query_start);
            let subject_span = 1 + record.subject_end.abs_diff(record.subject_start);
            let query_covered =
                (query_span as f64 / record.query_len as f64) >= MIN_ALIGN_COVER;
            let subject_covered =
                (subject_span as f64 / record.subject_len as f64) >= MIN_ALIGN_COVER;
            let equiv_len = min(query_span, subject_span) as f64
                >= (MIN_LEN_SIMILARITY * max(query_span, subject_span) as f64);

            (query_covered && subject_covered && equiv_len).then(|| {
                (
                    *seq_to_id.get(&record.query).unwrap(),
                    *seq_to_id.get(&record.target).unwrap(),
                )
            })
        })
        .collect();

    let graph = UnGraph::<(), (), usize>::from_edges(edges);
    let hash = _connected_components(&graph);
    hash.values()
        .map(|v| {
            v.iter()
                .flat_map(|id| id_to_seq.get(id).map(|v| v.to_string()))
                .collect::<Vec<_>>()
        })
        .collect()
}

// --------------------------------------------------
fn _connected_components<G>(graph: G) -> HashMap<usize, Vec<usize>>
where
    G: NodeCompactIndexable + IntoEdgeReferences + IntoNodeReferences,
{
    let mut vertex_sets = UnionFind::new(graph.node_bound());

    for edge in graph.edge_references() {
        let (a, b) = (edge.source(), edge.target());

        vertex_sets.union(graph.to_index(a), graph.to_index(b));
    }

    let mut components: HashMap<usize, Vec<usize>> = HashMap::new();

    for node in graph.node_references() {
        let node_idx = graph.to_index(node.id());
        let component_idx = vertex_sets.find(node_idx);

        components.entry(component_idx).or_default().push(node_idx);
    }

    components
}

// --------------------------------------------------
#[cfg(test)]
mod graph_tests {
    use crate::{graph::connected_components, parse_alignment};
    use anyhow::Result;
    use std::path::PathBuf;

    #[test]
    fn test_connected_components() -> Result<()> {
        let path = PathBuf::from("./tests/inputs/blast-consensi-self.tsv");
        let alignments = parse_alignment(&path)?;
        assert_eq!(alignments.len(), 100);

        let mut components = connected_components(alignments);
        components.sort_by_key(|v| v.len());

        // Two alignments were excluded by the filters
        assert_eq!(components.iter().map(|c| c.len()).sum::<usize>(), 98);

        // There are 47 groups
        assert_eq!(components.len(), 47);

        // The largest group has 46 members
        assert_eq!(components.last().unwrap().len(), 46);

        Ok(())
    }
}
