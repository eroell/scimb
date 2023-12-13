import pytest
import scanpy as sc

import scimb


def test_package_has_version():
    assert scimb.__version__ is not None


@pytest.mark.parametrize(
    "method",
    [
        "RandomOverSampler",
        #  "SMOTE"
    ],  # "SMOTENC", "SMOTEN", "ADASYN", "BorderlineSMOTE", "KMeansSMOTE", "SVMSMOTE"],
)
def test_methods_run(method):
    adata = sc.datasets.pbmc3k_processed()
    print(adata.obs.keys)

    scimb.pp.oversample(adata, key="louvain", method=method)


def test_method_not_implemented():
    adata = sc.datasets.pbmc3k_processed()
    print(adata.obs.keys)

    with pytest.raises(NotImplementedError):
        scimb.pp.oversample(adata, key="louvain", method="NotImplemented")


def test_data_type():
    adata = sc.datasets.pbmc3k_processed()
    print(adata.obs.keys)
    df = adata.obs.copy()

    with pytest.raises(ValueError):
        scimb.pp.oversample(df, key="louvain")
