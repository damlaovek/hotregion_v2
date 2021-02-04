from sklearn.cluster import DBSCAN
#from sklearn.cluster import OPTICS


class DoubleDBSCAN:
    def __init__(self, epsilon, min_points):
        self.epsilon = epsilon
        self.min_points = min_points
        
    def cluster(self, data):
        db = DBSCAN(eps=self.epsilon, min_samples=self.min_points).fit(data)
        #db = OPTICS(eps=self.epsilon, min_samples=self.min_points).fit(data)
        labels = db.labels_
        return labels
