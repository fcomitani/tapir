"""
...
F. Comitani     @2021
"""

def smart_selection(labs, to_select, how='any', val=1):

    """ Simplify selection of a single or multiple groups
        in a pandas dataframe.

        Args:
            labs (panda dataframe): one-hot-encoded classes membership dataframe 
                with samples as rows and classes as columns.
            to_select (list of int or strings): the groups to be selected.
            how (string): selection type, samples are chosen if they belong
                to 'any' or 'all' classes (default 'how').
            val (int): selection value, samples are chosen if their value in
                the classes memberhsip dataframe corresponds to this
                (default 1).
            
        Returns:
            (panda series): boolean series with the corresponding selection.
    """ 
    
    if isinstance(to_select,list):
        if how=='any':
            return (labs[to_select]==val).any(axis=1)
        elif how=='all':
            return (labs[to_select]==val).all(axis=1)
        else:
            return -1
    else:
        return labs[to_select]==val


def invert_dict(dictio):

    """ Invert dictionary whose values are lists and are not
        unique.

        Args:
            dictio (dictionary): the dictionary to invert.

        Returns:
            inv_dict (dictionary): a dictionary where the old values are 
                now keys and the old keys are now in the 
                value lists.
    """

    inv_dict={}

    for k, v in dictio.items():
        for g in v:
            if g not in inv_dict.keys():
                inv_dict[g] =   []
            inv_dict[g].append(k)

    return inv_dict