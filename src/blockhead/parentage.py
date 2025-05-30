def get_parentage(args):
    """
    Open the provided parentage file and get the F1, parent 1, and
    parent 2 sample names with p1/p2 in a tuple value per F1.
    """
    named_f1_dt = {}
    parent_dt = {}
    cross_ls = []

    with open(args.parentage) as f:
        for line in f:
            f1, p1, p2 = line.rstrip().split()
            named_f1_dt[f1] = (p1, p2)
            if p1 not in parent_dt:
                parent_dt[p1] = []
            if p2 not in parent_dt:
                parent_dt[p2] = []
            parent_dt[p1].append(f1)
            parent_dt[p2].append(f1)
            cross_ls.append((p1, p2))
    
    cross_ls = list(set(cross_ls))

    return parent_dt, cross_ls, named_f1_dt


def get_sample_ls(named_f1_dt: dict):
    """
    Iterate named_f1_dt and extract a unique list of samples.
    """
    sample_ls = []
    for k, v in named_f1_dt.items():
        sample_ls.append(k)
        for i in v:
            sample_ls.append(i)
    sample_ls = list(set(sample_ls))
    return sample_ls


def get_advanced(named_f1_dt):
    """
    Find the advanced hybrids based on parental tsv input.
    """
    adv_ls = []

    for k, v in named_f1_dt.items():
        for vi in v:
            if vi in named_f1_dt.keys():
                adv_ls.append(k)
    adv_ls = list(set(adv_ls))
    adv_ls.sort()
    return adv_ls


def get_advanced_lineage(adv, named_f1_dt):

    for i in named_f1_dt[adv]:
        if i in named_f1_dt.keys():
            f1 = i
            print(f"{i} is a hybrid")
            p1, p2 = named_f1_dt[f1]
            print(f"{p1} and {p2} are {f1} parents")
        else:
            p3 = i
            print(f"{i} is p3")
    
    lin_dt = {"adv" : adv,
              "f1"  : f1,
              "p1"  : p1,
              "p2"  : p2,
              "p3"  : p3}

    return lin_dt