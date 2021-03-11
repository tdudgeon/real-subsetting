# included fragments
match (m:V_XDSIP)-[:FRAG*]->(n:F2) where (n.chac > 0 or n.hac >= 5) and apoc.text.indexOf(n.smiles, ".") < 0 with distinct n as d return d.smiles;


# excluded fragments
match (m:V_XDSIP)-[:FRAG*]->(n:F2) where (n.chac = 0 and n.hac < 6) and apoc.text.indexOf(n.smiles, ".") < 0 with distinct n as d return d.smiles;

# included direct child fragments
match (m:V_XDSIP)-[:FRAG]->(n:F2) where ((n.chac > 0 and n.hac - n.chac > 1) or (n.chac = 0 and n.hac >= 5)) and apoc.text.indexOf(n.smiles, ".") < 0 with distinct n as d return d.smiles;


# dsip molecules
match (m:V_XDSIP) return m.smiles;



# running from shell
cat /tmp/query-direct.cql | cypher-shell -u neo4j -p <password> --format plain > /tmp/direct-children.smi
