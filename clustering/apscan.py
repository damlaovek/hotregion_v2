from sklearn.preprocessing import normalize
from sklearn.cluster import AffinityPropagation
import pandas as pd
import math
import numpy as np
from clustering.dbscan import DoubleDBSCAN

# Step 4: Merge Results
def calculate_argmax(result_i, j_indices, t):
    s = set(result_i)
    s.discard(-1)
    row = []
    for lbl_id in range(len(s)):
        counter = 0
        for index in j_indices:
            if result_i[index] == s[lbl_id]:
                counter += 1
        row.append[counter]
    arg_max = np.argmax(row)
    return s[arg_max]

def merge_two_clusters(result_i, result_j):
    merged_result = np.zeros(len(result_i))
    modified = np.zeros(len(result_i))
    merged_set = set(result_i) | set(result_j)
    new_cluster_label = max(merged_set) + 1
    for lbl in range(len(result_i)):
        if(modified[lbl] == 0):
            if result_i[lbl] == -1 and result_j[lbl] == -1:
                merged_result[lbl] = -1
                modified[lbl] = 1
            elif result_i[lbl] != -1 and result_j[lbl] != -1:
                merged_result[lbl] = result_i[lbl]
                modified[lbl] = 1
            elif result_i[lbl] == -1 and result_j[lbl] != -1:
                flag = False
                j_indices = []
                for i in range(len(result_j)):
                    if result_j[i] == result_j[lbl]:
                        j_indices.append(i)
                for i in range(len(result_j)):
                    if result_j[i] == result_j[lbl] and result_i[lbl] != -1:
                        arg_max = calculate_argmax(result_i, j_indices, result_j[lbl])
                        for j in j_indices:
                            merged_result[j] = arg_max
                            modified[j] = 1
                        flag = True
                if flag == False:
                    for j in j_indices:
                        merged_result[j] = new_cluster_label
                        modified[j] = 1
                    new_cluster_label += 1 
    return merged_result  


def merge_results(results):
    final_result = results[0]
    for result in range(1, len(results)):
        final_result = merge_two_clusters(final_result, results[result])
    return final_result

# Step 3: Double DBSCAN
def run_ddbscan(data, density_list):
    results = np.zeros((len(density_list)-1,data.shape[0]))
    for x in range(1, len(density_list)):
        min_points_new = density_list[0][0]/math.pow(density_list[x][1],2)*math.pow(density_list[0][1],2)
        eps2 = density_list[x][1]
        db = DoubleDBSCAN(eps2, min_points_new)
        results[x - 1] = db.cluster(data)
    return results
    
def run_ddbscan_version2(data, density_list):
    min_points_new = min(density_list, key = lambda t: t[0])[0]
    eps2 = sum(n for _, n, _ in density_list) / len(density_list)
    db = DoubleDBSCAN(eps2, min_points_new)
    results = db.cluster(data)
    return results

# Step 2: Normalized Density List Generation
def calc_distance(vec1, vec2):
	x1 = vec1[0]
	y1 = vec1[1]
	z1 = vec1[2]

	x2 = vec2[0]
	y2 = vec2[1]
	z2 = vec2[2]

	return math.sqrt(((x2-x1)**2)+((y2-y1)**2)+((z2-z1)**2))
	
def normalize_density_list(density_list):
    min_num = min(density_list, key = lambda t: t[0])[0]
    normalized_density_list = []
    for x in density_list:
        number_normalized = min_num
        radius_normalized = math.sqrt(number_normalized/x[0])*x[1]
        try:
            density_normalized = number_normalized/math.pow(radius_normalized, 2)
            normalized_density_list.append([number_normalized, radius_normalized, density_normalized])
        except OverflowError:
            pass
    return normalized_density_list

def generate_normalized_density_list(data, cluster_centers_indices, labels):
    n_clusters_ = len(cluster_centers_indices)
    density_list = []
    for k in range(n_clusters_):
        class_members = labels == k
        cluster_center = cluster_centers_indices[k]
        radius = 0
        x_center = data.loc[cluster_center]['X']
        y_center = data.loc[cluster_center]['Y']
        z_center = data.loc[cluster_center]['Z']
        distances = []
        num_members = 0
        for member in range(len(class_members)):
            if(class_members[member]):
                num_members += 1
                x_member = data.loc[member]['X']
                y_member = data.loc[member]['Y']
                z_member = data.loc[member]['Z']
                distance = calc_distance([x_center, y_center, z_center], [x_member, y_member, z_member])
                distances.append(distance)
                radius += distance
        radius = radius / num_members
        df = pd.DataFrame(distances, columns=['d'])
        df = df[df.d <= radius]
        if df.shape[0] > 1:
            number = df.shape[0]
            try:
                density = number/math.pow(radius, 2)
                density_list.append((number, radius, density))
            except OverflowError:
                pass
    density_list.sort(reverse=True, key=lambda tup: tup[2])
    #density_list = normalize_density_list(density_list)
    return density_list

# Step 1: Run Affinity Propagation    
def affinity_propagation(data):
    clusterer = AffinityPropagation()
    af = clusterer.fit(data)
    return af

def run_apscan(data):
    af = affinity_propagation(data)
    cluster_centers_indices = af.cluster_centers_indices_
    labels = af.labels_
    density_list = generate_normalized_density_list(data, cluster_centers_indices, labels)
    #cluster_results = run_ddbscan(data, density_list)
    #final_result = merge_results(cluster_results)
    final_result = run_ddbscan_version2(data, density_list)
    #print(set(final_result))
    return final_result
    

