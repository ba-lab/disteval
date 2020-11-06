import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import matplotlib.transforms as mtransforms
import matplotlib.text as mtext
from matplotlib.colors import LinearSegmentedColormap as lsCmap


# Don't currently have the ss generation code in so:
ss_dict_glbl = {}

# #######################################################################################
# #######################################################################################
# Complete Scatter Plot System by Jamie

## Numpy scheme helpers
def validate_dimensionality_and_return(ground_truth, prediction):
    gt_dim = ground_truth.ndim
    pd_dim = prediction.ndim

    # Assert dimensionality
    assert gt_dim == pd_dim, "ERROR: Ground truth and prediction dimensions do not match!"
    assert gt_dim < 3, "ERROR: Ground truth dimension >= 3 !"
    assert pd_dim < 3, "ERROR: Prediction dimension >= 3"
    assert gt_dim > 0, "ERROR: Ground truth dimension < 1 !"
    assert pd_dim > 0, "ERROR: Prediction dimension < 1 !"

    return ground_truth.ndim

def validate_sequence_length_and_return(ground_truth, prediction):
  assert ground_truth.shape == prediction.shape, "ERROR: Ground truth and prediction have different shapes!"
  return ground_truth.shape[0]

def get_flattened(dmap):
  if dmap.ndim == 1:
    return dmap
  elif dmap.ndim == 2:
    return dmap[np.triu_indices_from(dmap, k=1)]
  else:
    assert False, "ERROR: the passes array has dimension not equal to 2 or 1!"

def pre_process_for_metrics(ground_truth, prediction):
  _ = validate_sequence_length_and_return(ground_truth, prediction)
  _ = validate_dimensionality_and_return(ground_truth, prediction)

  gt_flat = get_flattened(ground_truth)
  pd_flat = get_flattened(prediction)
  return gt_flat, pd_flat

"""## Helpers"""

def get_separations(dmap):
  t_indices = np.triu_indices_from(dmap, k=1)
  separations = np.abs(t_indices[0] - t_indices[1])
  return separations
  
# return a 1D boolean array indicating where the sequence separation in the
# upper triangle meets the threshold comparison

def get_sep_thresh_b_indices(dmap, thresh, comparator):
  assert comparator in {'gt', 'lt', 'ge', 'le'}, "ERROR: Unknown comparator for thresholding!"
  dmap_flat = get_flattened(dmap)
  separations = get_separations(dmap)
  if comparator == 'gt':
    threshed = separations > thresh
  elif comparator == 'lt':
    threshed = separations < thresh
  elif comparator == 'ge':
    threshed = separations >= thresh
  elif comparator == 'le':
    threshed = separations <= thresh

  return threshed

# return a 1D boolean array indicating where the distance in the
# upper triangle meets the threshold comparison
def get_dist_thresh_b_indices(dmap, thresh, comparator):
  assert comparator in {'gt', 'lt', 'ge', 'le'}, "ERROR: Unknown comparator for thresholding!"
  dmap_flat = get_flattened(dmap)
  if comparator == 'gt':
    threshed = dmap_flat > thresh
  elif comparator == 'lt':
    threshed = dmap_flat < thresh
  elif comparator == 'ge':
    threshed = dmap_flat >= thresh
  elif comparator == 'le':
    threshed = dmap_flat <= thresh
  return threshed

"""## Function"""

# plot true distance (Y) against predicted distance (X) for locations where
def plot_distance_correlation_double_alpha(true_map, pred_map,
                                     min_sep=12, max_dist=[8, 20], 
                                     squeeze=True, color_sep=False,
                                     num_alpha_buckets=5, protein_name=None,
                                     pearson_8=0, pearson_20=0):
  s_thresh = min_sep
  d_thresh = max_dist
  # flattened
  true_map_flat = get_flattened(true_map)
  pred_map_flat = get_flattened(pred_map)

  # remove nans
  _nan_indices_true = np.isnan(true_map_flat)
  _nan_indices_pred = np.isnan(pred_map_flat)
  nan_indices = _nan_indices_true * _nan_indices_pred

  true_map_flat = true_map_flat[~nan_indices]
  pred_map_flat = pred_map_flat[~nan_indices]

  # boolean indexers
  b_indices_sep_gt_6 = get_sep_thresh_b_indices(true_map, s_thresh, 'gt')
    
  # separation colormap (if applicable)
  sep_cmap = lsCmap.from_list("",
                                [(0.00, "green"),
                                #  (0.05, "green"),
                                 (0.15, "xkcd:puke green"),
                                 (0.5, "xkcd:dark red"),
                                 (1.00, "xkcd:bright red")
                                 ]
                                )
  
  # get separations and remove nan
  separations = get_separations(true_map)
  separations = separations[~nan_indices]
  
  height = 8
  width = 8
  fig, ax = plt.subplots(1, 2, figsize=(width, height),
                         gridspec_kw={"width_ratios": d_thresh},
                         sharey=True)

  # setup plot 1
  # ##########################################################################
  b_indices_true_dist_8 = get_dist_thresh_b_indices(true_map, d_thresh[0], 'lt')

  # locations where the separation is > sep_thresh and the true distance < d_thresh
  b_indices_for_plot_1 = b_indices_sep_gt_6 & b_indices_true_dist_8

  # plot data
  true_data_1 = true_map_flat[b_indices_for_plot_1]
  pred_data_1 = pred_map_flat[b_indices_for_plot_1]

  err_data_1 = np.abs(true_data_1 - pred_data_1)

  # determine tick bounds
  true_max_1 = np.ceil(true_data_1.max()).astype(int)
  min_tick_1 = np.floor(min(true_data_1.min(), pred_data_1.min())).astype(int)
  max_tick_1 = np.ceil(max(true_data_1.max(), pred_data_1.max())).astype(int)

  # settings
  x_ticks_1 = range(min_tick_1, d_thresh[0]+1)

  # setup plot 2
  # ###################################################################################
  b_indices_true_dist_20 = get_dist_thresh_b_indices(true_map, d_thresh[1], 'lt')

  # locations where the separation is > sep_thresh and the true distance < dist_thresh
  b_indices_for_plot_2 = b_indices_sep_gt_6 & b_indices_true_dist_20

  # plot data
  true_data_2 = true_map_flat[b_indices_for_plot_2]
  pred_data_2 = pred_map_flat[b_indices_for_plot_2]

  err_data_2 = np.abs(true_data_2 - pred_data_2)

  err_bounds = np.linspace(-1, np.nanmax(err_data_2), num_alpha_buckets)

  # SETUP ALPHA
  alpha_list = np.geomspace(0.15, 10, num_alpha_buckets)


  # determine tick bounds
  true_max_2 = np.ceil(true_data_2.max()).astype(int)
  min_tick_2 = np.floor(min(true_data_2.min(), pred_data_2.min())).astype(int)
  max_tick_2 = np.ceil(max(true_data_2.max(), pred_data_2.max())).astype(int)

  # settings
  # x_ticks_2 = range(min_tick_2, d_thresh[1]+1)
  
  # PLOT FOR  D < 8A
  # ###################################################################################
  if color_sep:
    separations_1 = separations[b_indices_for_plot_1]

    sc_1 = ax[0].scatter(true_data_1, pred_data_1, 
                    alpha=0.4, c=err_data_1, cmap=sep_cmap, vmax=20)
  else:
    ax[0].scatter(true_data_1, pred_data_1, alpha=0.55)

  ax[0].plot([min_tick_1, true_max_1+1], [min_tick_1, true_max_1+1], ls='--', c='k')
  ax[0].set_xlabel("True distance", fontsize=14)
  ax[0].set_ylabel("Predicted distance", fontsize=14)
  # _ = ax[0].set_xticks(x_ticks_1)
  # _ = ax[0].set_title(f"True < {d_thresh[0]} : r = {get_pearson(true_data_1, pred_data_1):.2f}",
  #                  fontsize=16)
  _ = ax[0].text(s=f"r = {pearson_8:.4f}",
                   fontsize=14, x=4, y=2)
  
  #######################


  # PLOT FOR  D < 20A
  # ###################################################################################
  if color_sep:
    separations_2 = separations[b_indices_for_plot_2]
    for i in range(len(err_bounds) - 1):
      _plot_data_indices = (err_data_2 > err_bounds[i]) & (err_data_2 <= err_bounds[i+1])
      _true_plot_data = true_data_2[_plot_data_indices]
      _pred_plot_data = pred_data_2[_plot_data_indices]
      _err_plot_data = err_data_2[_plot_data_indices]

      sc_2 = ax[1].scatter(_true_plot_data, _pred_plot_data,
                        alpha=alpha_list[i], c=_err_plot_data, cmap=sep_cmap, vmin=0, vmax=20
                        )

    # sc_2 = ax[1].scatter(true_data_2, pred_data_2, 
    #                 alpha=0.4, c=err_data_2, cmap=sep_cmap, vmax=20)
  else:
    ax[1].scatter(true_data_2, pred_data_2, alpha=0.55)
    



  ax[1].plot([min_tick_2, true_max_2+1], [min_tick_2, true_max_2+1], ls='--', c='k')
  ax[1].set_xlabel("True distance", fontsize=14)
  # _ = ax[1].set_xticks(x_ticks_2)
  # _ = ax[1].set_title(f"True < {d_thresh[1]} : r = {get_pearson(true_data_2, pred_data_2):.2f}",
  #                  fontsize=16)
  _ = ax[1].text(s=f"r = {pearson_20:.4f}",
                   fontsize=14, x=14, y=2)
  
  if color_sep:
    # colorbar is same for both
    cb = fig.colorbar(sc_2,
                        ticks=np.arange(2, 22, 2).astype(int),
                        label="Absolute Error (in Å)", 
    )

  # ########## GLOBAL ########

  # ############################################
  # y ticks are shared between the two
  y_ticks = range(1, 30)

  _ = ax[0].set_yticks(y_ticks)
  _ = ax[1].set_yticks(y_ticks)

  title_str = "Residue separations ≥ 12"
  if protein_name is not None:
    title_str += f" - {protein_name}"
  fig.suptitle(title_str,
               y=1.015, fontsize=16)

  fig.tight_layout()
  return fig

# #######################################################################################
# #######################################################################################
# #### Complete heatmap system by Jamie

### Generic helpers

def clip_distance_map(in_map, lower_thresh, upper_thresh):
  _map = in_map.copy()
  _map[in_map > upper_thresh] = upper_thresh
  _map[in_map < lower_thresh] = lower_thresh

  return _map

# just check if the argument is valid
# put in tiny wrapper since it is used often
def check_upper_lower(_arg):
  assert _arg in {"upper", "lower"}

# Determine if a map is binary, that is, the unique non-nan values are 0, 1
def is_binary(in_map):
  unq = np.unique(in_map)
  unq = unq[~np.isnan(unq)]
  unq = set(unq)
  return unq == {0, 1}

# Get a binary contact map from a generic map
# Accepts: R-valued distances, contact probabilities, and binary contact maps
def get_binary_contact(in_map, thresh=None):
  map_bin = in_map.copy()

  if np.nanmax(in_map) > 2:
    if thresh is None:
      _thr = 8
    else:
      _thr = Thresh
    map_bin[in_map < _thr] = 1
    map_bin[in_map >= _thr] = 0

  else:
    if not is_binary(map_bin):
      if thresh is None:
        _thr = 0.5
      else:
        _thr = Thresh
      map_bin[in_map >= _thr] = 1
      map_bin[in_map < _thr] = 0
      
  return map_bin

# Return maps of confusion metrics in this order:
#   true_positive_map, false_positive_map, true_negative_map, false_negative_map
def get_confusion_maps(true_map, pred_map):
  true_binary = get_binary_contact(true_map)
  pred_binary = get_binary_contact(pred_map)

  tp_map = np.full_like(true_binary, np.nan)
  fp_map = np.full_like(true_binary, np.nan)

  tn_map = np.full_like(true_binary, np.nan)
  fn_map = np.full_like(true_binary, np.nan)

  tp_map[(pred_binary == 1) & (true_binary == 1)] = 1.
  fp_map[(pred_binary == 1) & (true_binary == 0)] = 1.

  tn_map[(pred_binary == 0) & (true_binary == 0)] = 1.
  fn_map[(pred_binary == 0) & (true_binary == 1)] = 1

  return tp_map, fp_map, tn_map, fn_map

def get_absolute_error_map(true_map, pred_map):
  return np.abs(true_map - pred_map)

def get_log10_error_map(true_map, pred_map):
  _err = get_absolute_error_map(true_map, pred_map)
  _err = np.log10(_err + 1)
  _err[_err < 0] = 0
  return _err

# blank the triangle NOT in the used_triangle argument
def blank_unused_triangle(_map, used_triangle):
  check_upper_lower(used_triangle)
  indices = np.tril_indices_from(_map)
  _out_map = _map.copy()
  _out_map[indices] = np.nan
  
  if used_triangle == "lower":
    _out_map = _out_map.T

  return _out_map

a = np.arange(5).astype(float)
b = np.arange(3, 8).astype(float)
a[2] = np.nan
a-b

"""### Build helpers"""

# build semantic color maps from list of form
# [
#   (np.nan, mp_color, bool=label this Amstrong tick),
#   (int: Amstrongs, mpl_color, bool=label this Amstrong tick),
# ...,
# (int: Amstrongs, mpl_color, bool=label this Amstrong tick),
# (np.inf: Amstrongs, mpl_color, bool=label this Amstrong tick)
# ]
# Note: np.nan and np.inf are special values for map minimum and maximum

def build_cmap_from_param_list(_cmap_param_list, _map):
  norm_cols = lambda l, h, x: (x - l) / (h - l)
  
  map_min = np.nanmin(_map)
  map_max = np.nanmax(_map)
  
  col_list = [(num, col) for num, col, _ in _cmap_param_list[1:-1] if col is not None]
  col_list.insert(0, (map_min, _cmap_param_list[0][1]))
  col_list.append((map_max, _cmap_param_list[-1][1]))
  col_list = [(norm_cols(map_min, map_max, c), n) for (c, n) in col_list]

  dist_cmap = lsCmap.from_list("", col_list)
  return dist_cmap

# Wrappper to render a map in a single triangle of an axis
def render_map(_map, _cmap, _axis, _vmin, _vmax, upper_or_lower):
  check_upper_lower(upper_or_lower)
  _cm = matplotlib.cm.get_cmap(_cmap)

  if _vmin is None:
    _min = np.nanmin(_map)
  else:
    _min = _vmin
  if _vmax is None:
    _max = np.nanmax(_map)
  else:
    _max = _vmax

  map_to_render = blank_unused_triangle(_map, used_triangle=upper_or_lower)
  _ims = _axis.imshow(map_to_render, _cm, vmin=_min, vmax=_max,
                      aspect='auto')#, interpolation="antialiased")

  return _ims

# call as one of final steps
# must construct ticks on the map axis first
def render_orthogonal_markers(_L, _axis, stride, upper_or_lower):
  check_upper_lower(upper_or_lower)

  # place at every other tick
  if upper_or_lower == "upper":
    _ticks = _axis.get_xticks()
  else:
    _ticks = _axis.get_yticks()
  _orth_ticks = list(range(0, _L, stride))
  _orth_ticks.pop(0)

  # line style
  _line_params = {
      "ls": '-',
      "color": "blue",
      "lw": 2,
      "alpha": 0.3
  }

  # render
  for t in _orth_ticks:
    if upper_or_lower == "upper":
      _axis.axhline(t, t/_L, 1, **_line_params)
      _axis.axvline(t, (1-t/_L), 1, **_line_params)
    else:
      _axis.axhline(t, 0, t/_L, **_line_params)
      _axis.axvline(t, 0, (1-t/_L), **_line_params)


# diagonal lines at 6, 12, 24
def render_diagonal_markers(_L, _axis, upper_or_lower):
  check_upper_lower(upper_or_lower)

  # line style
  _line_params = {
      "ls": '--',
      "lw": 2,
      "color": "blue",
      "alpha": 0.3
  }
  

  # render
  for t in [6, 12, 24]:
    if upper_or_lower == "upper":
      _x_bounds = [t, _L-1]
      _y_bounds = [1, _L - t]
    else:
      _x_bounds = [1, _L - t]
      _y_bounds = [t, _L-1]
    _axis.plot(_x_bounds, _y_bounds, **_line_params)

def build_right_color_bar(_up_ims, _axis, _ticks, _tick_lbls):
  _vert_cbar = plt.colorbar(_up_ims, cax=_axis, 
                            ticks=_ticks)
  _vert_cbar.ax.set_yticklabels(_tick_lbls)

  return _vert_cbar

def build_bottom_color_bar(_dn_ims, _axis, _ticks, _tick_lbls):
  _hori_cbar = plt.colorbar(_dn_ims, cax=_axis, 
                            ticks=_ticks, orientation="horizontal")
  _hori_cbar.ax.set_xticklabels(_tick_lbls)
  _hori_cbar.ax.invert_xaxis()

  return _hori_cbar

def build_sequence_color_bar(_L, _axis, _cmap, _width, _ticks, _prot_name):
  _cm = matplotlib.cm.get_cmap(_cmap)
  _seq_arr = np.arange(0, _L)
  _seq_arr = np.broadcast_to(_seq_arr, (_width, _L))

  _seq_cbar = _axis.imshow(_seq_arr, _cm, aspect="auto")

  _seq_title = "Color for residue indices"
  if _prot_name is not None:
    _seq_title += f" of {_prot_name}"

  _axis.set_title(_seq_title, fontsize=14)
  _axis.set_xticks(_ticks)
  _axis.tick_params(axis='y', which='both', left=False, labelleft=False)
  _axis.tick_params(axis='x', which='both', bottom=True, labelbottom=False)

  return _seq_cbar

# put the sequence index color map on the main diagonal.  
# without this main diagonal is black.  may cover first couple neighbor diagonals
def render_main_diagonal(_L, _axis, _ss_dict, _cmap):
  _diag_map = matplotlib.cm.get_cmap(_cmap)
  for i in range(_L-1):
    _color = _diag_map(float(i) / _L)
    if i != 0:
      if i+1 in _ss_dict and _ss_dict[i+1] == 'H': _color = 'red'
      elif i+1 in _ss_dict and _ss_dict[i+1] == 'E': _color = 'green'
      elif i+1 in _ss_dict and _ss_dict[i+1] == 'C': _color = "#eeeeee"
    _axis.plot([i, i+1], [i, i+1],
              #  linewidth = max([-0.025 * _L + 8.5, 8]), 
               linewidth=8,
               alpha = 1, color=_color
               )

# various ticks and labels for various purposes
def get_strided_ticks(_L, _stride):
  _ticks = list(range(0, _L + 1, _stride))
  _tick_lbls = [str(t) for t in _ticks]
  return _ticks, _tick_lbls

def get_probability_ticks():
  _ticks = [0.0, 0.25, 0.5, 0.75, 1.0]
  _tick_lbls = [f"{t:.2f}" for t in _ticks]
  return _ticks, _tick_lbls

def get_absolute_error_ticks():
  _ticks = [0, 5, 10, 15, 20]
  _tick_lbls = [str(t) for t in _ticks]
  return _ticks, _tick_lbls

def get_log10_error_ticks():
  _ticks = np.arange(0, 1.25, 0.1)
  _tick_lbls = [f"{t:.2f}" for t in _ticks]
  return _ticks, _tick_lbls

def get_standard_distance_ticks(_map):
  _ticks = [np.nanmin(_map).min(), 6, 8, 10, 12, 14, 16, 18, 20]
  _tick_labels = [f"{_ticks[0]:.2f}", "6", "8", '', "12", '', "16", '', "20+"]
  return _ticks, _tick_labels

# select the ticks
def find_ticks(kind=None, _map=None):
  if kind is None:
    return None
  elif kind in {"distance", "true_distance", "pred_distance"}:
    return get_standard_distance_ticks(_map)
  elif kind in {"probability", "pred_contact"}:
    return get_probability_ticks()
  elif kind == "absolute_error":
    return get_absolute_error_ticks()
  elif kind == "log10_error":
    return get_log10_error_ticks()
  else:
    return None

# select the label for x/y of main heatmap axis
def find_map_label(kind=None, sub_kind=None):
  if kind == None:
    return None
  elif kind == "distance":
    return "Distance (in Å)"
  elif kind == "true_distance":
    return "True Distance (in Å)"
  elif kind == "pred_distance":
    return "Predicted Distance (in Å)"
  elif kind == "true_contact":
    return "True Contacts"
  elif kind in {"probability", "pred_contact"}:
    return "Predicted Contact Probabilities"
  elif kind == "absolute_error":
    return "Absolute Error (in Å)"
  elif kind == "log10_error":
    return r"log$_{10}($Absolute Error$)$"
  elif kind == "confusion":
    return "Confusion Map\n(green = TP, red = FP, blue = FN"
  else:
    return None

# construct ticks and labels from a semantic color list as above
def get_cbar_ticks_and_labels(_cmap_param_list, _map, capped=True):
  map_min = np.nanmin(_map)
  map_max = np.nanmax(_map)
  
  # get ticks
  _ticks = [n for n, _, _ in _cmap_param_list]
  _ticks[0] = map_min
  _ticks[-1] = map_max
  _ticks = [map_min, 6, 8, 10, 12, 14, 16, 18, map_max]

  # get labels
  _tick_lbls = [str(n) if b else '' for n, _, b in _cmap_param_list]
  _tick_lbls[0] = f"{map_min:.3f}"
  max_lbl = str(int(map_max))
  if capped:
    max_lbl += '+'
  _tick_lbls[-1] = max_lbl

  return _ticks, _tick_lbls

# put ticks and labels on the heatmap, color main diagonal black
def decorate_heatmap(_L, _axis, _ticks, _tick_lbls, 
                     _right_label, _bottom_label,
                     use_upper=True, use_lower=True):
  _axis.set_xticks(_ticks)
  _axis.set_yticks(_ticks)

  _identity = np.eye(_L)
  _identity[_identity == 0] = np.nan

  _axis.imshow(_identity, "binary_r", aspect="auto")

  _axis.tick_params(which='both',
                    left=use_lower, labelleft=use_lower,
                    right=False, labelright=False,
                    top=use_upper, labeltop=use_upper,
                    bottom=False, labelbottom=False
                    )
  
  if _right_label is not None:
    _axis.yaxis.set_label_position("right")
    _axis.set_ylabel(_right_label, labelpad=16)
  
  if _bottom_label is not None:
    _axis.set_xlabel(_bottom_label, labelpad=16)

  if not use_upper:
    _axis.spines["top"].set_visible(False)
    _axis.spines["right"].set_visible(False)

  if not use_lower:
    _axis.spines["bottom"].set_visible(False)
    _axis.spines["left"].set_visible(False)



"""### Inner Heatmap Build Function"""

def build_full_heatmap(upper_map=None, lower_map=None,
                       upper_type=None, lower_type=None,
                       protein_name=None,
                       seq_cmap=None,
                       lower_cmap=None, upper_cmap=None,
                       u_vmin=None, u_vmax=None,
                       l_vmin=None, l_vmax=None,
                       use_right_bar=False,
                       use_bottom_bar=False,
                       use_main_diagonal=False,
                       ss_dict=None
                       ):

  # derive sequence length
  if upper_map is not None:
    seq_len = len(upper_map)
  if lower_map is not None:
    seq_len = len(lower_map)


  # Get upper map color map
  _uvmin = u_vmin
  _uvmax = u_vmax
  if upper_cmap is not None:
    # case for the semantic coloring list
    if isinstance(upper_cmap, list):
      u_cm = build_cmap_from_param_list(upper_cmap, upper_map)
      _uvmin = None
      _uvmax = None
    else:
      u_cm = matplotlib.cm.get_cmap(upper_cmap)
      _uvmin = u_vmin
      _uvmax = u_vmax
    upper_label = find_map_label(upper_type)

  # Get lower map color map
  _lvmin = l_vmin
  _lvmax = l_vmax
  if lower_cmap is not None:
    # case for the semantic coloring list
    if isinstance(lower_cmap, list):
      l_cm = build_cmap_from_param_list(lower_cmap, lower_map)
      _lvmin = None
      _lvmax = None
    else:
      l_cm = matplotlib.cm.get_cmap(lower_cmap)

    lower_label = find_map_label(lower_type)


  # ############################################################################

  # Build Figure and Gridspec
  hm_fig = plt.figure(figsize=(8, 8), constrained_layout=True)
  _wr = None
  _hr = [2, 96]
  _nbars = 0
  _ncols = 1
  _nrows = 2
  if use_right_bar:
    _wr = [96, 4]
    _ncols += 1
  if use_bottom_bar:
    _hr += [4]
    _nrows += 1
  hm_gs = hm_fig.add_gridspec(_nrows, _ncols,
                              width_ratios=_wr,
                              height_ratios=_hr,
                              wspace=0, hspace=0
                              )
  
  # Place sequence and heatmap subplots
  sbar_ax = hm_fig.add_subplot(hm_gs[0, 0])
  cmap_ax = hm_fig.add_subplot(hm_gs[1, 0], sharex=sbar_ax)

  # Calculate tick strides from length
  if seq_len < 100:
    tick_stride = 10
  elif seq_len <= 200:
    tick_stride = 20
  elif seq_len <= 400:
    tick_stride = 40
  elif seq_len <= 1000:
    tick_stride = 100
  elif seq_len <= 2000:
    tick_stride = 200
  elif seq_len <= 5000:
    tick_stride = 500
  else:
    tick_stride = 1000

  # Get ticks for sequence / heatmap
  idx_ticks, idx_tick_lbls = get_strided_ticks(seq_len, tick_stride)

  # Decorate heatmap - ticks, labels
  if upper_type is not None:
    _use_upp = True
  else:
    _use_upp = False
  if lower_type is not None:
    _use_low = True
  else:
    _use_low = False

  r_lbl = find_map_label(upper_type)
  b_lbl = find_map_label(lower_type)
  decorate_heatmap(seq_len, cmap_ax, idx_ticks, idx_tick_lbls,
                   _right_label=r_lbl, _bottom_label=b_lbl,
                   use_upper=_use_upp, use_lower=_use_low
                   )

  # Place additional components
  build_sequence_color_bar(seq_len, sbar_ax, seq_cmap, 2, idx_ticks, protein_name)
  if upper_type is not None:
    if upper_type == "absolute_error":
      _up_map = get_absolute_error_map(upper_map, lower_map)
    elif upper_type == "log10_error":
      _up_map = get_log10_error_map(upper_map, lower_map)
    elif upper_type == "confusion":
      _up_map = upper_map
      _tp, _fp, _, _fn = get_confusion_maps(lower_map, upper_map)
      _ = render_map(_tp, "summer_r", cmap_ax, 0, 1, "upper")
      _ = render_map(_fp, "autumn_r", cmap_ax, 0, 1, "upper")
      _ = render_map(_fn, "cool_r", cmap_ax, 0, 1, "upper")
    else:
      _up_map = upper_map
    if upper_type != "confusion":
      upper_ims = render_map(_up_map, u_cm, cmap_ax, _uvmin, _uvmax, "upper")

    # Markers
    render_orthogonal_markers(seq_len, cmap_ax, tick_stride*2, "upper")
    render_diagonal_markers(seq_len, cmap_ax, "upper")
  else:
    lower_ims = None

  # Lower heatmap
  if lower_type is not None:
    if lower_type == "absolute_error":
      _low_map = get_absolute_error_map(upper_map, lower_map)
    elif lower_type == "log10_error":
      _low_map = get_log10_error_map(upper_map, lower_map)
    elif lower_type == "confusion":
      _low_map = lower_map
      _tp, _fp, _, _fn = get_confusion_maps(lower_map, upper_map)
      _ = render_map(_tp, "summer_r", cmap_ax, 0, 1, "lower")
      _ = render_map(_fp, "autumn_r", cmap_ax, 0, 1, "lower")
      _ = render_map(_fn, "cool_r", cmap_ax, 0, 1, "lower")
    else:
      _low_map = lower_map
    if lower_type != "confusion":
      lower_ims = render_map(_low_map, l_cm, cmap_ax, _lvmin, _lvmax, "lower")

    # Markers
    render_orthogonal_markers(seq_len, cmap_ax, tick_stride*2, "lower")
    render_diagonal_markers(seq_len, cmap_ax, "lower")
  else:
    lower_ims = None

  # Handle rightside bar
  if use_right_bar:
    rbar_ax = hm_fig.add_subplot(hm_gs[1, 1])
    right_ticks, right_tick_lbls = find_ticks(upper_type, upper_map)
    _ = build_right_color_bar(upper_ims, rbar_ax,
                              right_ticks, right_tick_lbls
                              )
  # Handle bottom bar
  if use_bottom_bar:
    bbar_ax = hm_fig.add_subplot(hm_gs[2, 0])
    bottom_ticks, bottom_tick_lbls = find_ticks(lower_type, lower_map)
    _ = build_bottom_color_bar(lower_ims, bbar_ax,
                               bottom_ticks, bottom_tick_lbls
                               )
    
  if use_main_diagonal:
    render_main_diagonal(seq_len, cmap_ax, ss_dict, seq_cmap)
  
  return hm_fig

"""## Build parameter packs and find plot
* use this to determine your parameters from two simple arguments
* then unpack the parameter dictionary in the build heatmap function
* the protein_name is optional, if set to None it won't be printed
"""

# PARAMETER PACKS
def get_parameter_pack(plot_type, color_scheme='standard', protein_name=None):
  if plot_type is None:
    return None

  true_dmap_cols = [
                  (np.nan, 'xkcd:light khaki', True),
                  (6, 'xkcd:light khaki', True),
                  (8, 'xkcd:dusty orange', True),
                  (10, 'xkcd:easter purple', False),
                  (12, "xkcd:robin's egg", True),
                  (14, 'xkcd:yellowish green', False),
                  (16, None, True),
                  (18, None, False),
                  (np.inf, 'xkcd:cool grey', True)
                  ]
                  
  pred_dmap_cols = [
                    (np.nan, 'xkcd:light khaki', True),
                    (6, 'xkcd:light khaki', True),
                    (8, 'xkcd:faded orange', True),
                    (10, 'xkcd:pale violet', False),
                    (12, 'xkcd:pale aqua', True),
                    (14, 'xkcd:pale olive green', False),
                    (16, None, True),
                    (18, None, False),
                    (np.inf, 'xkcd:cool grey', True)
                    ]

  # Three color schemes
  if color_scheme == "standard":
    true_d_cm = "YlGn_r"
    pred_d_cm = "OrRd_r"
    prob_c_cm = "Reds"
    log_er_cm = "gray_r"
    abs_er_cm = "Reds"

    seq_id_cm = "jet"

  elif color_scheme == "semantic":
    true_d_cm = true_dmap_cols
    pred_d_cm = pred_dmap_cols
    prob_c_cm = "Reds"
    log_er_cm = "gray_r"
    abs_er_cm = "Reds"

    seq_id_cm = "jet"

  elif color_scheme == "uniform":
    true_d_cm = "plasma_r"
    pred_d_cm = "viridis_r"
    prob_c_cm = "cividis_r"
    log_er_cm = "cividis_r"
    abs_er_cm = "inferno_r"

    seq_id_cm = "jet"



  # Distance vs. Distance
  _pp_dist_vs_dist = dict(
    upper_type="pred_distance", lower_type="true_distance",
    seq_cmap=seq_id_cm,
    lower_cmap=true_d_cm, upper_cmap=pred_d_cm,
    u_vmin=None, u_vmax=20,
    l_vmin=None, l_vmax=20,
    use_right_bar=True,
    use_bottom_bar=True,
    use_main_diagonal=True,
    protein_name=protein_name
  )

  # Contact probabilities vs. true distance
  _pp_con_vs_prob = dict(
    upper_type="pred_contact", lower_type="true_distance",
    seq_cmap=seq_id_cm,
    lower_cmap=true_d_cm, upper_cmap=prob_c_cm,
    u_vmin=0.0, u_vmax=1.0,
    l_vmin=None, l_vmax=20,
    use_right_bar=True,
    use_bottom_bar=True,
    use_main_diagonal=True,
    protein_name=protein_name
  )

  # Only true distance
  _pp_true_dist = dict(
    upper_type="true_distance", lower_type=None,
    seq_cmap=seq_id_cm,
    lower_cmap=None, upper_cmap=true_d_cm,
    u_vmin=None, u_vmax=20,
    l_vmin=None, l_vmax=None,
    use_right_bar=True,
    use_bottom_bar=False,
    use_main_diagonal=True, 
    protein_name=protein_name 
  )

  # Only predicted distance
  _pp_pred_dist = dict(
    upper_type="pred_distance", lower_type=None,
    seq_cmap=seq_id_cm,
    lower_cmap=None, upper_cmap=pred_d_cm,
    u_vmin=None, u_vmax=20,
    l_vmin=None, l_vmax=None,
    use_right_bar=True,
    use_bottom_bar=False,
    use_main_diagonal=True,  
    protein_name=protein_name
  )

  # Only contact probabilities
  _pp_pred_prob = dict(
    upper_type="probability", lower_type=None,
    seq_cmap=seq_id_cm,
    lower_cmap=None, upper_cmap=prob_c_cm,
    u_vmin=0, u_vmax=1,
    l_vmin=None, l_vmax=None,
    use_right_bar=True,
    use_bottom_bar=False,
    use_main_diagonal=True,  
    protein_name=protein_name
  )

  # Log error vs. confusion
  _pp_err_vs_conf = dict(
      upper_type="log10_error", lower_type="confusion",
      seq_cmap=seq_id_cm,
      lower_cmap=None, upper_cmap=log_er_cm,
      u_vmin=0, u_vmax=1.25,
      l_vmin=0, l_vmax=1,
      use_right_bar=True,
      use_bottom_bar=False,
      use_main_diagonal=True,
      protein_name=protein_name
  )

  # Log error vs. absolute error
  _pp_log_vs_abs = dict(
      upper_type="log10_error", lower_type="absolute_error",
      seq_cmap=seq_id_cm,
      lower_cmap=abs_er_cm, upper_cmap=log_er_cm,
      u_vmin=0, u_vmax=1.25,
      l_vmin=0, l_vmax=None,
      use_right_bar=True,
      use_bottom_bar=True,
      use_main_diagonal=True,
      protein_name=protein_name
  )

  # Log error vs. true distance
  _pp_err_vs_dist = dict(
      upper_type="log10_error", lower_type="distance",
      seq_cmap=seq_id_cm,
      lower_cmap=true_d_cm, upper_cmap=log_er_cm,
      u_vmin=0, u_vmax=1.25,
      l_vmin=None, l_vmax=20,
      use_right_bar=True,
      use_bottom_bar=True,
      use_main_diagonal=True,
      protein_name=protein_name
  )

  # Only Log error
  _pp_log_error = dict(
      upper_type="log10_error", lower_type=None,
      seq_cmap=seq_id_cm,
      lower_cmap=None, upper_cmap=log_er_cm,
      u_vmin=0, u_vmax=1.25,
      l_vmin=None, l_vmax=None,
      use_right_bar=True,
      use_bottom_bar=False,
      use_main_diagonal=True,
      protein_name=protein_name
  )

  # Only Confusion
  _pp_confusion = dict(
      upper_type="confusion", lower_type=None,
      seq_cmap=seq_id_cm,
      lower_cmap=None, upper_cmap=None,
      u_vmin=0, u_vmax=1,
      l_vmin=None, l_vmax=None,
      use_right_bar=False,
      use_bottom_bar=False,
      use_main_diagonal=True,
      protein_name=protein_name
  )

  if plot_type == "distance_vs_distance":
    pp_pack = _pp_dist_vs_dist
  elif plot_type in {"distance_vs_probability", "distance_vs_contact"}:
    pp_pack = _pp_con_vs_prob
  elif plot_type == "only_pred_distance":
    pp_pack = _pp_pred_dist
  elif plot_type == "only_true_distance":
    pp_pack = _pp_true_dist
  elif plot_type == "only_probability":
    pp_pack = _pp_pred_prob
  elif plot_type == "log_error_vs_confusion":
    pp_pack = _pp_err_vs_conf
  elif plot_type == "log_error_vs_abs_error":
    pp_pack = _pp_log_vs_abs
  elif plot_type == "only_log_error":
    pp_pack = _pp_log_error
  elif plot_type == "only_confusion":
    pp_pack = _pp_confusion
  elif plot_type == "log_error_vs_true_distance":
    pp_pack = _pp_err_vs_dist

  return pp_pack


# #######################################################################################
# #######################################################################################
# ########### CHORD DIAGRAMS
def dmap2chordimage(ND = None, chord_file = None, ss={}, cmap_for_chord='jet'):
    chord_fig, chord_ax = plt.subplots(1, 1, figsize=(16, 16))
    cb_map = np.copy(ND)
    L = len(ND)
    size_of_each_protein_in_chord_map = 360 / L
    rnum2angle = {}
    rnum2xy = {}
    rnum2indexxy = {}
    num_2_index_xy = {}
    for i in range(L + 1):
        a1 = round((i-1) * 2 * math.pi / L, 5)
        a2 = round(i * 2 * math.pi / L, 5)
        x = round(0.50 + 0.495 * math.cos((a1+a2)/2), 5)
        y = round(0.50 + 0.50 * math.sin((a1+a2)/2), 5)
        rnum2angle[str(i) + "a1"] = str(a1)
        rnum2angle[str(i) + "a2"] = str(a2)
        rnum2xy[str(i) + "x"] = x
        rnum2xy[str(i) + "y"] = y
        lx = round(0.47 + 0.52 *math.cos((a1+a2)/2) + 0.04 * math.cos(i * 2 * math.pi / L), 5)
        ly = round(0.49 + 0.52 * math.sin((a1+a2)/2) + 0.02 * math.sin(i * 2 * math.pi / L), 5)
        rnum2indexxy[str(i) + "x"] = lx
        rnum2indexxy[str(i) + "y"] = ly
        num_2_index_xy[i-1] = {"x": lx, "y": ly}
    cb_map[np.isnan(cb_map)] = 1000.0
    cb_map[cb_map < 4.0] = 4.0
    cb_map = 4.0 / cb_map


    ## ######### ADDED BY JAMIE

    cmap = matplotlib.cm.get_cmap(cmap_for_chord)

    ## ######### DONE WITH ADDITION BY JAMIE

    pair_n_intensity = {}
    for i in range(L):
        for j in range(i, L):
            # Only show the strong connections
            if cb_map[i, j] < 0.1: continue
            # Skip local connections
            if abs(i - j ) < 6: continue
            pair_n_intensity[str(i+1) + ' ' + str(j+1)] = cb_map[i, j] * cb_map[i, j]
    tokeep = {}
    keepcount = 0
    for pair, intensity in reversed(sorted(pair_n_intensity.items(), key = lambda x: x[1])):
        tokeep[pair + ' ' + str(intensity)] = int(pair.split()[0])
        keepcount += 1
        if keepcount > 5 * L: break
    for pair_n_intensity, ii in sorted(tokeep.items(), key = lambda x: x[1]):
        i, j, intensity = pair_n_intensity.split()
        intensity = float(intensity)

        chord_ax.plot([rnum2xy[str(i) + "x"], rnum2xy[str(j) + "x"]], [rnum2xy[str(i) + "y"], rnum2xy[str(j) + "y"]], 
          linewidth = 6 * intensity, color= cmap(float(i) / L), alpha = intensity, zorder = -1)

    for i in range(L):
        mycolor = cmap(float(i) / L)
        if i == 0: mycolor = cmap(float(i) / L)

        e1 = matplotlib.patches.Arc(xy=(.5, .5), width=1, height=1, linewidth=12, angle=i * size_of_each_protein_in_chord_map, 
          theta2=size_of_each_protein_in_chord_map, color=mycolor)
        chord_ax.add_patch(e1)

        step = int(L/20)
        if i % step != 0: continue
        if i > L - step/2: continue
        chord_ax.text(num_2_index_xy[i]["x"], num_2_index_xy[i]["y"], i+1, fontsize=32, color= cmap(float(i) / L))

    for i in range(L):
        mycolor = cmap(float(i) / L)
        if i+1 in ss and ss[i+1] == 'H': mycolor = 'red'
        if i+1 in ss and ss[i+1] == 'E': mycolor = 'green'
        if i+1 in ss and ss[i+1] == 'C': mycolor = "#eeeeee"
        if i == 0: mycolor = cmap(float(i) / L)
        e1 = matplotlib.patches.Arc(xy=(.5, .5), width=1, height=1, linewidth=7, angle=i * size_of_each_protein_in_chord_map, 
          theta2=size_of_each_protein_in_chord_map, color=mycolor)

        chord_ax.add_patch(e1)

        step = int(L/20)
        if i % step != 0:
            continue
        if i > L - step/2:
            continue
        chord_ax.text(num_2_index_xy[i]["x"], num_2_index_xy[i]["y"], i+1, fontsize=32, color= cmap(float(i) / L))
    chord_ax.axis('off')
    return chord_fig


# #######################################################################################
# #######################################################################################
# Create the visualizations

def do_visualization(script_args, native_basename, basename, ND, D, NC, C, pearson_8=0, pearson_20=0):
  # Cases:    Have ND & D:
#               heatmap:    distance_vs_distance
#               errormap:   log_error_vs_true_distance
#           Have ND & C:
#               heatmap:    distance_vs_probability
#               errormap:   only_confusion
#           Have only D (or ND):
#               heatmap:    only_pred_distance (or only_true_distance)
#               errormap:   None
#           Have C:
#               heatmap:    only_pred_prob
#               errormap:   None
#           Have only NC:   *not possible / supported*

    # native_basename : filename for native map
    # basename : filename for prediction
    # ss : filename for secondary structure

    print("\n\n")

    if ND is not None:
      ND_clipped = clip_distance_map(ND, 3.5, 20)
    else:
      ND_clipped = None
    if D is not None:
      D_clipped = clip_distance_map(D, 3.5, 20)
    else:
      D_clipped = None

    # Decode colorscheme argument
    if   script_args.color_scheme == 1: _col_scheme = "standard"
    elif script_args.color_scheme == 2: _col_scheme = "semantic"
    elif script_args.color_scheme == 3: _col_scheme = "uniform"

    # Heatmap parameter packs
    hm_pack = None
    er_pack = None

    # Filenames
    hm_file_name = None
    er_file_name = None

    if ND_clipped is not None:
        if D_clipped is not None:
            hm_low_map = ND_clipped
            hm_upp_map = D_clipped

            hm_pack = get_parameter_pack(plot_type="distance_vs_distance", color_scheme=_col_scheme)
            er_pack = get_parameter_pack(plot_type="log_error_vs_abs_error", color_scheme=_col_scheme)

            hm_file_name = native_basename + '.vs.' + basename + '.heatmap.' + script_args.vis_filetype
            er_file_name = native_basename + '.vs.' + basename + '.errormap.' + script_args.vis_filetype

            print("Visualizing true distances vs. predicted distances, saving to:", hm_file_name)
            print("Visualizing errors as log(absolute_error) vs. absolute error, saving to:", er_file_name)
        elif C is not None:
            hm_low_map = ND_clipped
            hm_upp_map = C

            hm_pack = get_parameter_pack(plot_type="distance_vs_probability", color_scheme=_col_scheme)
            er_pack = get_parameter_pack(plot_type="only_confusion", color_scheme=_col_scheme)

            hm_file_name = native_basename + '.vs.' + basename + '.heatmap.' + script_args.vis_filetype
            er_file_name = native_basename + '.vs.' + basename + '.errormap.' + script_args.vis_filetype
            
            print("Visualizing true distances vs. predicted contacts, saving to:", hm_file_name)
            print("Visualizing errors as confusion map, saving to:", er_file_name)
        else:
            hm_upp_map = ND_clipped
            hm_low_map = None

            hm_pack = get_parameter_pack(plot_type="only_true_distance", color_scheme=_col_scheme)
            
            hm_file_name = native_basename + '.vs.' + basename + '.heatmap.' + script_args.vis_filetype

            print("Visualizing true distances, saving to:", hm_file_name)
    elif D_clipped is not None:
        hm_upp_map = D_clipped
        hm_low_map = None

        hm_pack = get_parameter_pack(plot_type="only_pred_distance", color_scheme=_col_scheme)


        hm_file_name = native_basename + '.vs.' + basename + '.heatmap.' + script_args.vis_filetype

        print("Visualizing predicted distances, saving to:", hm_file_name)
    elif C is not None:
        hm_upp_map = C
        hm_low_map = None

        hm_pack = get_parameter_pack(plot_type="only_pred_prob", color_scheme=_col_scheme)

        hm_file_name = native_basename + '.vs.' + basename + '.heatmap.' + script_args.vis_filetype

        print("Visualizing predicted contacts, saving to:", hm_file_name)
    else:
        print("Nothing to visualize.")

    # Make heatmap plots
    if hm_pack is not None:
        hm_fig = build_full_heatmap(upper_map=hm_upp_map, lower_map = hm_low_map, ss_dict=ss_dict_glbl, **hm_pack)
        hm_fig.savefig(hm_file_name, bbox_inches='tight', format=script_args.vis_filetype)

    if er_pack is not None:
        er_fig = build_full_heatmap(upper_map=hm_upp_map, lower_map = hm_low_map, ss_dict=ss_dict_glbl, **er_pack)
        er_fig.savefig(er_file_name, bbox_inches='tight', format=script_args.vis_filetype)


    # Make chord diagrams
    c_seq_cmap = hm_pack["seq_cmap"]
    if ND is not None:

      native_chord_fig_name = native_basename + '.chord_diagram.' + script_args.vis_filetype 
      native_chord_fig = dmap2chordimage(ND=ND, ss=ss_dict_glbl, cmap_for_chord=c_seq_cmap)

      native_chord_fig.savefig(native_chord_fig_name, bbox_inches='tight', format=script_args.vis_filetype)

      print("Saving chord diagram: predicted, to:", native_chord_fig_name)


    if D is not None:
      
      pred_chord_fig_name = basename + '.chord_diagram.' + script_args.vis_filetype
      pred_chord_fig = dmap2chordimage(ND=D, ss=ss_dict_glbl, cmap_for_chord=c_seq_cmap)

      pred_chord_fig.savefig(pred_chord_fig_name, bbox_inches='tight', format=script_args.vis_filetype)

      print("Saving chord diagram: predicted, to:", pred_chord_fig_name)

    # Make scatter plots
    if (ND is not None) and (D is not None):
      
      scatter_fig_name = native_basename + '.vs.' + basename + '.scatterplot.' + script_args.vis_filetype
      scatter_fig = plot_distance_correlation_double_alpha(ND, D,
                                     min_sep=12, max_dist=[8, 20], 
                                     squeeze=True, color_sep=True,
                                     num_alpha_buckets=5, protein_name=None,
                                     pearson_8=pearson_8, pearson_20=pearson_20)
      scatter_fig.savefig(scatter_fig_name, bbox_inches='tight', format=script_args.vis_filetype)

      print("Saving scatterplot to:", scatter_fig_name)