import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import numpy as np

def validate_clustering(clustering: list, num_vertices):
  assert num_vertices > 3
  if (len(clustering)==0):
    return
  alreadySeen = np.zeros((num_vertices,), dtype=bool)
  for cluster_idx, members in enumerate(clustering):
      assert len(set(members)) == len(members), f"cluster ID {cluster_idx} contains duplicate city"
      for city_id in members:
          if city_id >= num_vertices or city_id < 0:
              raise ValueError(f"Invalid city ID {city_id} in cluster {cluster_idx}")
          if alreadySeen[city_id] == 1:
              raise ValueError(f"City ID {city_id} already appeared in another cluster (of a lower index)")
          alreadySeen[city_id] = 1
  # technically we don't need this any more at this point but just to be clear
  assert np.all(alreadySeen==1), "some cities are not assigned to any clusters"

def evaluate_vertex2cluster(valid_clustering: list, num_vertices: int) -> np.ndarray:
  assert len(valid_clustering) >= 2
  v2c = - np.ones((num_vertices,), dtype=int)
  for cluster_id, members in enumerate(valid_clustering):
    for v in members:
      v2c[v] = cluster_id
  return cluster_id

def validate_permutation(seq: np.ndarray, num_vertices: int, elem_name = "vertex"):
  seq = seq.squeeze()
  assert seq.ndim == 1
  assert len(seq) == num_vertices
  assert np.all(seq >= 0)
  assert np.all(seq < num_vertices)
  assert len(set(seq)) == len(seq), f"some {elem_name} is visited more than once"
  

def extract_cluster_seq_and_validate(valid_clustering: list, tour: np.ndarray, num_vertices: int):
  assert np.all(tour >= 0)
  assert np.all(tour < num_vertices)
  if len(valid_clustering) == 0:
    return [1]
  v2c = evaluate_vertex2cluster(valid_clustering, num_vertices)
  cluster_seq = []
  for vertex_id in tour:
    cluster_seq.append(v2c[vertex_id])
  validate_permutation(cluster_seq, len(valid_clustering), "cluster")
  return cluster_seq
    

def eval_cost(
  graph_xyz: np.ndarray, tour: np.ndarray, 
  clustering: list, rounding = True, ord = 2) -> float:
  """
  works for both geometric TSP and geometric generalized TSP and 
  """
  assert graph_xyz.ndim == 2
  validate_clustering(clustering, len(graph_xyz))
  is_generalized = len(clustering) != 0
  expected_tour_len = len(clustering) if is_generalized else len(graph_xyz) 
  assert len(tour) == expected_tour_len, f"|tour| = {len(tour)}, expect = {expected_tour_len}"
  
  tourXyz = graph_xyz[[*tour, tour[0]]]
  edgeCosts = np.linalg.norm(np.diff(tourXyz, axis=0), axis=1, ord=ord)
  if rounding: # the TSPLIB convention
    edgeCosts = np.round(edgeCosts, decimals=0)
  return np.sum(edgeCosts)
  
  

def plot_tour(ax: Axes, graph_xy: np.ndarray, tour: np.ndarray):
  tour = tour.flatten()
  assert len(tour) == len(graph_xy)
  ax.plot(*graph_xy[[*tour, tour[0]]].T, '-or', ms=3)
  ax.set_aspect('equal')
  
if __name__ == "__main__":
  import sys
  import pathlib
  import logging
  logging.getLogger().setLevel(logging.INFO)
  thisDir = pathlib.Path(__file__).parent
  sys.path.append(thisDir)
  from tsplib import read_tsplib_tour, read_tsplib_geomTsp
  
  from argparse import ArgumentParser
  p = ArgumentParser("Validate and evaluate a TSP tour")
  p.add_argument("prob_file", type=str)
  p.add_argument("tour_file", type=str)
  args = p.parse_args()
  
  logging.info(f"Reading problem file: {args.prob_file}")
  xy, clustering, prob_name = read_tsplib_geomTsp(args.prob_file)
  is_generalized = len(clustering) != 0
  logging.info(f" NAME: {prob_name}")
  logging.info(f" node coordinates: {xy.shape[0]}x{xy.shape[1]}")
  if is_generalized:
    logging.info(f"number of clusters: {len(clustering)}")
    validate_clustering(clustering)
  
  logging.info(f"Reading tour file: {args.tour_file}")
  tour, tour_name = read_tsplib_tour(args.tour_file)
  logging.info(f" NAME: {tour_name}")
  if is_generalized:
    cluster_seq = extract_cluster_seq_and_validate(clustering, tour, len(xy))
    logging.info(f"Cluster sequence: {cluster_seq}")
  else:
    validate_permutation(tour, len(xy))
  
  cost = eval_cost(xy, tour, clustering)
  logging.info(f"Tour cost {cost:.3f}")
  
  if xy.shape[1] != 2:
    logging.warning("Skipping visualization because it supports only 2D TSP, not {xy.shape[1]}-D")
  else:
    _, ax = plt.subplots()
    if is_generalized:
      # first, let's plot the clusters
      for members in clustering:
        ax.plot(xy[members],"-xb")
    plot_tour(ax, xy, tour)
    ax.set_title(f"{prob_name} (cost: {cost:.0f})")
    plt.show()
  