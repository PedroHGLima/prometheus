import numpy as np 
from scipy.special import kl_div

def rebin_histogam(old_bin_edges, old_bin_count, new_bin_edges):
    '''
    This function is used to rebin a numpy histogram given new edges.

    Arguments:

    - old_bin_edges: the old edged to be replaced;
    - old_bin_counts: the old counts to be reorganized;
    - new_bin_edges: the new edged to be used.
    '''
    x1 = np.asarray(old_bin_edges)
    y1 = np.asarray(old_bin_count)
    x2 = np.asarray(new_bin_edges)

    # the fractional bin locations of the new bins in the old bins
    i_place = np.interp(x2, x1, np.arange(len(x1)))

    cum_sum = np.r_[[0], np.cumsum(y1)]

    # calculate bins where lower and upper bin edges span
    # greater than or equal to one original bin.
    # This is the contribution from the 'intact' bins (not including the
    # fractional start and end parts.
    whole_bins = np.floor(i_place[1:]) - np.ceil(i_place[:-1]) >= 1.
    start = cum_sum[np.ceil(i_place[:-1]).astype(int)]
    finish = cum_sum[np.floor(i_place[1:]).astype(int)]

    y2 = np.where(whole_bins, finish - start, 0.)

    bin_loc = np.clip(np.floor(i_place).astype(int), 0, len(y1) - 1)

    # fractional contribution for bins where the new bin edges are in the same
    # original bin.
    same_cell = np.floor(i_place[1:]) == np.floor(i_place[:-1])
    frac = i_place[1:] - i_place[:-1]
    contrib = (frac * y1[bin_loc[:-1]])
    y2 += np.where(same_cell, contrib, 0.)

    # fractional contribution for bins where the left and right bin edges are in
    # different original bins.
    different_cell = np.floor(i_place[1:]) > np.floor(i_place[:-1])
    frac_left = np.ceil(i_place[:-1]) - i_place[:-1]
    contrib = (frac_left * y1[bin_loc[:-1]])

    frac_right = i_place[1:] - np.floor(i_place[1:])
    contrib += (frac_right * y1[bin_loc[1:]])

    y2 += np.where(different_cell, contrib, 0.)

    return y2

def combine_edges( counts, edges, threshold ):
    '''
    This function will merge the bins into a numpy histogram given a threshold

    Arguments:

    - counts: the histogram bin counts.
    - edges: the histogram bin edges.
    - threshold: the minimum counts for each bin.
    '''

    max_ix = counts.argmax()
    c_list = list( counts )   # Lists can be popped from
    e_list = list( edges )    # Lists can be popped from

    def eliminate_left( ix ):
        # Sum the count and eliminate the edge relevant to ix
        # Before the peak (max_ix)
        nonlocal max_ix
        max_ix -= 1         # max_ix will change too.
        c_list[ix+1]+=c_list[ix]
        c_list.pop(ix)
        e_list.pop(ix+1)

    def eliminate_right( ix ):
        # Sum the count and eliminate the edge relevant to ix
        # after the peak (max_ix) 
        c_list[ix-1]+=c_list[ix]
        c_list.pop(ix)
        e_list.pop(ix)

    def first_lt():
        # Find the first ix less than the threshold
        for ix, ct in enumerate( c_list[:max_ix] ):
            if ct < threshold:
                return ix  # if ct < threshold return the index and exit the function
        # The function only reaches here if no ct values are less than the threshold
        return -1  # If zero items < threshold return -1

    def last_lt():
        # Find the last ix less than the threshold
        for ix, ct in zip( range(len(c_list)-1, max_ix, -1), c_list[::-1]):
            # ix reduces from len(c_list)-1, c_list is accessed in reverse order.
            if ct < threshold:
                return ix
        return -1  # If no items < threshold return -1

    cont = True
    while cont:
        # Each iteration removes any counts less than threshold
        # before the peak.  This process would combine e.g. counts of [...,1,2,3,...] into [..., 6, ...]
        ix = first_lt()
        if ix < 0:
            cont = False   # If first_lt returns -1 stop while loop
        else:
            eliminate_left( ix )

    cont = True
    while cont:
        ix = last_lt()
        if ix < 0:
            cont = False   # If last_lt returns -1 stop while loop
        else:
            eliminate_right( ix )

    return np.array( c_list ), np.array( e_list )

def calc_chi_s(counts_obs_hist, counts_ref_hist):
    '''
    This function will calculate the $\chi^s$ statistic for a given two histograms

    Arguments:

    - counts_obs_hist: coints of observed numpy histogram.
    - counts_ref_hist: counts of reference numpy histogram.
    '''

    return (counts_obs_hist - counts_ref_hist)/np.sqrt(counts_ref_hist)

def calc_chi_s_residuals(counts_obs_hist, counts_ref_hist):
    '''
    This function will calculate the $\chi^s$ residuals for $\chi^s$ statistic.

    Arguments:

    - counts_obs_hist: coints of observed numpy histogram.
    - counts_ref_hist: counts of reference numpy histogram.

    '''
    diff_counts = (counts_obs_hist.sum() - counts_ref_hist.sum())
    return (diff_counts*(counts_ref_hist/counts_ref_hist.sum()))/np.sqrt(counts_ref_hist)

def calc_kl(pk, qk):
    '''
    A function to calculate the Kullback-Libler divergence between p and q distribution.
    Arguments:
    pk: pdf values from p distribution
    qk: pdf values from q distribution
    '''
    return kl_div(pk, qk)

def calc_js(pk, qk):
    '''
    A function to calculate the Jensen-Shanon divergence between p and q distribution.
    Arguments:
    pk: pdf values from p distribution
    qk: pdf values from q distribution
    '''
    mk = 0.5*(pk+qk)
    return 0.5*(calc_kl(pk, mk) + calc_kl(qk, mk))