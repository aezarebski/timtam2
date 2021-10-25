graph TD
    origin[Origin]
    style origin fill:#deebf7
    tree[Tree]
    style tree fill: #deebf7
    params[RealParameter*]
    style params fill:#deebf7
    twpp[TreeWithPointProcess]
    style twpp fill:#3182bd
    timtam[TimTam]
    style timtam fill:#3182bd
    points[TraitSet]
    style points fill:#deebf7

    origin --> twpp
    origin --> timtam
    twpp --> timtam
    points --> twpp
    params --> timtam
    tree --> twpp
    tree --> timtam
