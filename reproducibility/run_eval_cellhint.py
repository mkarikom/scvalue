import scanpy as sc
import cellhint

methods = ['Uniform' 'GeoSketch' 'Sphetcher' 'Hopper' 'KH' 'scSampler' 'scValue']

cols = [
    'Cycling T&NK',
    'Tnaive/CM_CD4_activated',
    'Trm_gut_CD8',
    'MAIT',
    'Trm/em_CD8',
    'Tnaive/CM_CD8',
    'Tem/emra_CD8',
    "Tgd_CRTAM+",
    'Tregs',
    'Trm_Tgd',
]

sketch_size = 20066 # 10% of all datasets
for method in methods:
    adata = sc.read('experiments/cellhint_spleen/%s.%d.h5ad' % (method, sketch_size))
    print(adata.obs.Dataset.value_counts())

    alignment = cellhint.harmonize(adata, 'Dataset', 'Original_annotation')
    print(alignment)

    path = 'experiments/cellhint_spleen/spleen_alignment.%s.%d.pkl' % (method, sketch_size)
    alignment.write(path)
    alignment = cellhint.DistanceAlignment.load(path)
    print(alignment.relation.head(10))

    df = alignment.relation.copy()
    col_names = df.columns
    if method in ['sphetcher', 'kh']:
        df = df[df[col_names[4]].isin(cols)].copy()
    else:
        df = df[df[col_names[0]].isin(cols)].copy()

    cellhint.treeplot(df, save='experiments/cellhint_spleen/%s_T.%d.png' % (method, sketch_size))
