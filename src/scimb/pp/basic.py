from typing import Union

from anndata import AnnData
from imblearn.over_sampling import (
    ADASYN,
    SMOTE,
    SMOTEN,
    SMOTENC,
    SVMSMOTE,
    BorderlineSMOTE,
    KMeansSMOTE,
    RandomOverSampler,
)


def basic_preproc(adata: AnnData) -> int:
    """Run a basic preprocessing on the AnnData object.

    Parameters
    ----------
    adata
        The AnnData object to preprocess.

    Returns
    -------
    Some integer value.
    """
    print("Implement a preprocessing function here.")
    return 0


# TODO: Add an oversampling function which calls the imblearn oversampling library
# and can do oversampling, taking as arguments
# - an AnnData object
# - a key of anndata's .obs as the class to balance
# - the oversampling name (e.g. SMOTE)
# - the oversampling parameters (e.g. k_neighbors=5) as **kwargs
# It should return a new AnnData object.


def oversampling(
    adata: AnnData, key: str, method: str = "RandomOverSampler", use_rep: str = None, **kwargs
) -> Union[AnnData, None]:
    """_summary_

    Args:
        adata (AnnData): _description_
        key (str): _description_
        method (str): _description_
        use_rep (str, optional): _description_. Defaults to None.

    Raises
    ------
        ValueError: _description_
        NotImplementedError: _description_

    Returns
    -------
        Union[AnnData, None]: _description_
    """
    # TODO: check input object

    data_is_adata = isinstance(adata, AnnData)

    if not data_is_adata:
        raise ValueError(f"Input data is not an AnnData object: type of {adata}, is {type(adata)}")

    # ensure that random state is set to 0 if not set by user
    if kwargs is None:
        kwargs = {"random_state": 0}
    elif "random_state" not in kwargs.keys():
        kwargs["random_state"] = 0

    if method == "RandomOverSampler":
        sampler = RandomOverSampler(**kwargs)
    elif method == "SMOTE":
        sampler = SMOTE(**kwargs)
    elif method == "SMOTENC":
        sampler = SMOTENC(**kwargs)
    elif method == "SMOTEN":
        sampler = SMOTEN(**kwargs)
    elif method == "ADASYN":
        sampler = ADASYN(**kwargs)
    elif method == "BorderlineSMOTE":
        sampler = BorderlineSMOTE(**kwargs)
    elif method == "KMeansSMOTE":
        sampler = KMeansSMOTE(**kwargs)
    elif method == "SVMSMOTE":
        sampler = SVMSMOTE(**kwargs)
    else:
        raise ValueError(f"Unknown oversampling method: {method}")

    if use_rep is None:
        use_data = adata.X
    elif use_rep in adata.obsm.keys():
        use_data = adata.obsm[use_rep]
    else:
        raise ValueError(f"Error with use_rep: is not None and is not in adata.obsm: {use_rep}.")

    if key in adata.obs.keys():
        use_label = adata.obs[key]
    else:
        raise ValueError(f"Error with key: is not in adata.obs: {key}.")

    # TODO: only interested in index, that is sampler.sample_indices_
    # is then fit_resample really the way to go?
    # _, y = sampler.fit_resample(use_data, use_label)
    sampler.fit_resample(use_data, use_label)

    # get the index of the new data
    # idx = y.index # not useful it just gets a range index

    # sample the adata
    # adata_sub = adata[idx, :].copy()
    adata_sub = adata[sampler.sample_indices_, :].copy()

    # store sampling information in adata.obs
    adata_sub.obs["orig_pos"] = sampler.sample_indices_
    adata_sub.obs["orig_index"] = adata.obs.index[sampler.sample_indices_]

    # TODO: think if we want to keep the original data. Representations might become misleading
    # First thought: keep .obs, .var, .uns (e.g. to keep the log1p base info) (but not rank genes etc?)
    # drop .obsm, .varm, .obsp, any tracks of pca and neighbors and umaps
    # TODO: might keep some things - user responsibility & decision how to proceed?
    # I personally would probably only keep the count data and the .obs and .var but depends on the use case?
    # TODO: think whether to rename duplicated indices (only an upsampling thing)

    return adata_sub


# TODO: Add an undersampling function which calls the imblearn undersampling library
# and can do undersampling, taking as arguments
# - an AnnData object
# - a key of anndata's .obs as the class to balance
# - the oversampling name (e.g. SMOTE)
# - the oversampling parameters (e.g. k_neighbors=5) as **kwargs
# It should return a new AnnData object.
# Potentially this could be merged with the oversampling function above.
def undersample(adata: AnnData, key: str, method: str, **kwargs) -> Union[AnnData, None]:
    """_summary_

    Args:
        adata (AnnData): _description_
        key (str): _description_
        method (str): _description_

    Raises
    ------
        NotImplementedError: _description_

    Returns
    -------
        AnnData: _description_
    """
    raise NotImplementedError("Implement undersampling here.")


# TODO: Add a function which calls the imblearn combination library
# y
