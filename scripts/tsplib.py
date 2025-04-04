from enum import Enum
import numpy as np

class _Problem_parse_state(Enum):
    seeking_name = 0
    seeking_type = 1
    seeking_num_cities = 2
    seeking_num_clusters = 3
    seeking_edge_type = 4
    seeking_NODE_COORD_SECTION = 5
    seeking_GTSP_SET_SECTION = 6
    
class _Tour_parse_state(Enum):
    seeking_name = 0
    seeking_type = 1
    seeking_dim = 2
    seeking_tour = 3

def read_tsplib_geomTsp(filepath: str) -> tuple[np.ndarray, list[list[int]], str]:
    """
    @return the xy data
    @return membership (only relevant for generalized TSP)
    @return name of the problem
    """
    state = _Problem_parse_state.seeking_name
    is_generalized = False
    with open(filepath, 'r') as f:
        line_original = f.readline() 
        line = line_original.upper().strip()
        while line:
            if state == _Problem_parse_state.seeking_name:
                if line.startswith("NAME"):
                    if not ":" in line:
                        raise ValueError("missing colon `:` in the row 'NAME' ")
                    name =  line_original.split(":")[1].strip()
                    state = _Problem_parse_state.seeking_type
            elif state == _Problem_parse_state.seeking_type:
                if line.startswith("TYPE"): 
                  if line.endswith("TSP\n"):
                    state = _Problem_parse_state.seeking_num_cities
                  elif line.endswith("GTSP\n"):
                    state = _Problem_parse_state.seeking_num_cities
                    is_generalized = True
                  else:
                    raise ValueError("Unrecognized "+line)         
            elif state == _Problem_parse_state.seeking_num_cities:
                if line.startswith("DIMENSION"):
                    num_cities = int(line.split(':')[-1].strip())
                    state = _Problem_parse_state.seeking_num_clusters if \
                        is_generalized else \
                          _Problem_parse_state.seeking_edge_type
            elif state == _Problem_parse_state.seeking_num_clusters:
                if line.startswith("GTSP_SETS"):
                    num_clusters = int(line.split(':')[-1].strip())
                    state = _Problem_parse_state.seeking_edge_type
            elif state == _Problem_parse_state.seeking_edge_type:
                if line.startswith("EDGE_WEIGHT_TYPE"):
                    weight_type = line.split(':')[-1].strip()
                    # print("Weight type:", weight_type)
                    try:
                        num_geom_dim = int(weight_type[-2])
                    except ValueError:
                        raise ValueError(
                            f"Expecting a geometric TSP/GTSP problem but got EDGE_WEIGHT_TYPE: {weight_type}")
                    state = _Problem_parse_state.seeking_NODE_COORD_SECTION
            elif state == _Problem_parse_state.seeking_NODE_COORD_SECTION:
                if line.startswith("NODE_COORD_SECTION"):
                    xy = np.zeros((num_cities, num_geom_dim), dtype=float)
                    for ii in range(num_cities):
                        d = list(map(float, f.readline().split())) # get rid of whitespaces and the last character "\n"
                        if len(d) != 1 + num_geom_dim:
                            raise ValueError(f"Invalid NODE_COORD_SECTION: row {ii+1} should have {num_geom_dim+1} items, got {len(d)} instead")
                        if d[0] != ii+1:
                            raise ValueError(f"Invalid NODE_COORD_SECTION: row {ii+1} is mis-numbered as {d[0]}")
                        xy[ii] = d[1:]
                    if is_generalized:                    
                      state = _Problem_parse_state.seeking_GTSP_SET_SECTION
                    else:
                      return xy, [], name
            elif state == _Problem_parse_state.seeking_GTSP_SET_SECTION:
                if line.startswith("GTSP_SET_SECTION"):
                    clusters = []
                    for ii in range(num_clusters):
                        d = list(map(int, f.readline().split())) # get rid of whitespaces and the last character "\n"
                        if d[-1] != -1:
                            raise ValueError(f"Invalid GTSP_SET_SECTION: row {ii+1} missing the -1 delimiter")
                        if d[0] != ii+1:
                            raise ValueError(f"Invalid GTSP_SET_SECTION: row {ii+1} is mis-numbered as {d[0]}")
                        # members of this cluster (note: -1 because we use zero-indexing)
                        members = np.array(d[1:-1],dtype=int)-1
                        clusters.append(members) 
                    return xy, clusters, name
            line = f.readline()
    
    raise ValueError("Invalid GTSPLIB file because it fails to progress while " + str(state.name))

def read_tsplib_tour(filepath) -> tuple[np.ndarray, str]:
    """parse a tour file and convert it to 0-based indexing

    Expected format of the file:
    ```
    NAME : myTour
    COMMENT : blablabla
    TYPE : TOUR
    DIMENSION : 34523 (length of the tour)
    TOUR_SECTION
    ..
    ..
    -1
    EOF
    ```
    
    Note that (1) no preceding space in each line and 
    (2) the order matters    
    
    @return: tour
    @return: name
    """   
    state = _Tour_parse_state.seeking_name
    tour_length = 0
    with open(filepath, 'r') as f:
        line_original = f.readline() 
        line = line_original.upper()
        while line:
            if state == _Tour_parse_state.seeking_name:
                if line.startswith("NAME"):
                    if not ":" in line:
                        raise ValueError("missing colon `:` in the row 'NAME' ")
                    name =  line_original.split(":")[1].strip()
                    state = _Tour_parse_state.seeking_type
            elif state == _Tour_parse_state.seeking_type:
                if line.startswith("TYPE") and line.endswith("TOUR\n"):
                    state = _Tour_parse_state.seeking_dim
            elif state == _Tour_parse_state.seeking_dim:
                if line.startswith("DIMENSION"):
                    tour_length = int(line.split(' ')[-1])
                    state = _Tour_parse_state.seeking_tour
            elif state == _Tour_parse_state.seeking_tour:
                if line.startswith("TOUR_SECTION"):
                    tour_data = np.zeros((tour_length,), dtype=int)
                    for ii in range(tour_length):
                        d = int(f.readline()[:-1]) # get rid of the last character "\n"
                        if d <= 0:
                            raise ValueError(f"Invalid TOUR_SECTION: containing erroneous city index {d}")
                        tour_data[ii] = d
                    # reinstate 0-indexing
                    tour_data -= 1
                    
                    # penultimate line
                    pline = f.readline()
                    if pline != "-1\n":
                        raise ValueError("Invalid TOUR_SECTION: longer than the declared dimension")
                    
                    if len(tour_data) != len(set(tour_data)):
                        raise ValueError("Invalid TOUR_SECTION: some city is visited more than once.")
                    
                    return tour_data.reshape(-1), name
                
            line = f.readline()
    
    raise ValueError("Invalid tour file because it fails to progress while " + str(state.name))