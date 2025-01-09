use std::collections::HashMap;

use petgraph::{
    graph::UnGraph,
    unionfind::UnionFind,
    visit::{EdgeRef, IntoEdgeReferences, IntoNodeReferences, NodeCompactIndexable, NodeRef},
};

use crate::RmBlastOutput;

/// Given a slice of `RmBlastOutput` records, get the connected components of a graph where edges
/// are created between sequences that have been aligned.
///
/// The components are returned as `Vec<Vec<&str>>`, where each inner `Vec` contains the names of the
/// sequences in the component.
pub fn connected_components<'a>(records: &'a [RmBlastOutput]) -> Vec<Vec<&'a str>> {
    let mut seq_to_id = HashMap::new();
    let mut id_to_seq = HashMap::new();

    let mut id_cnt = 0usize;

    records.iter().for_each(|record| {
        for &s in &[record.query, record.target] {
            seq_to_id.entry(s).or_insert_with(|| {
                let id = id_cnt;
                id_to_seq.insert(id, s);
                id_cnt += 1;
                id
            });
        }
    });

    let edges: Vec<_> = records
        .iter()
        .map(|r| {
            (
                *seq_to_id.get(r.query).unwrap(),
                *seq_to_id.get(r.target).unwrap(),
            )
        })
        .collect();

    let graph = UnGraph::<(), (), usize>::from_edges(edges);
    let components = _connected_components(&graph);

    components
        .values()
        .map(|v| {
            v.iter()
                .map(|id| *id_to_seq.get(id).unwrap())
                .collect::<Vec<&str>>()
        })
        .collect()
}

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
