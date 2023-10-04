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


def oversample(adata: AnnData, key: str, method: str, **kwargs) -> Union[AnnData, None]:
    """_summary_

    Args:
        adata (AnnData): _description_
        key (str): _description_
        method (str): _description_

    Raises
    ------
        ValueError: _description_
        NotImplementedError: _description_

    Returns
    -------
        Union[AnnData, None]: _description_
    """
    if method == "RandomOverSampler":
        RandomOverSampler(**kwargs)
    if method == "SMOTE":
        SMOTE(**kwargs)
    elif method == "SMOTENC":
        SMOTENC(**kwargs)
    elif method == "SMOTEN":
        SMOTEN(**kwargs)
    elif method == "ADASYN":
        ADASYN(**kwargs)
    elif method == "BorderlineSMOTE":
        BorderlineSMOTE(**kwargs)
    elif method == "KMeansSMOTE":
        KMeansSMOTE(**kwargs)
    elif method == "SVMSMOTE":
        SVMSMOTE(**kwargs)
    else:
        raise ValueError(f"Unknown oversampling method: {method}")

    raise NotImplementedError("Implement oversampling here.")

    return


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
