import json

sys.path.append('/Users/cathalye/Packages/histoannot/')
import phas.client.api as phas
from phas.dltrain import spatial_transform_roi


def connect_to_server(task):
    # PHAS connection parameters
    PHAS_URL='https://histo.itksnap.org'
    PRIVATE_KEY = '/Users/cathalye/.private/histoitk_api_key.json'

    conn = phas.Client(PHAS_URL, PRIVATE_KEY)
    # Create a sampling ROI task object to pass to Slide class for downloading sampling ROI json
    task = phas.SamplingROITask(conn, task)

    return task


def read_json_property(jsonfile, property):
    with open(jsonfile, 'r') as file:
        data = json.load(file)
        slide_id = data.get(f'{property}', None)
    return slide_id


def get_nearest_chunk_map(chunk_mask):
    # Explanation of this function in nearest_chunk_map.ipynb
    chunk_mask_arr = sitk.GetArrayFromImage(chunk_mask)

    chunk_labels = np.unique(chunk_mask_arr)
    chunk_labels = chunk_labels[chunk_labels != 0] # Remove the background label

    def _get_dist_from_chunk(k):
        # Get the distance of every pixel from the boundary of the chunk with label k
        mask = sitk.BinaryThreshold(chunk_mask, int(k), int(k), 1, 0) # Extract only the chunk with label k
        return sitk.SignedDanielssonDistanceMap(mask, insideIsPositive=False, squaredDistance=True)


    dist_maps_all_chunks = { k: _get_dist_from_chunk(k) for k in chunk_ids }
    dist_map_all_chunks = np.array([ dist_maps_all_chunks[k] for k in chunk_ids ])
    # Find the chunk with the minimum distance to a given point
    nearest_chunk = chunk_ids[np.argmin(dist_maps_all_chunks, axis = 0)]

    nearest_chunk_map = np.reshape(nearest_chunk, chunk_mask_arr.shape)

    return nearest_chunk_map
