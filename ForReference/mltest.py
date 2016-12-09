import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error

file_path = 'D:/Project/python_learning/games.csv'
# Read in the data.
games = pd.read_csv(file_path)
# Print the names of the columns in games.
print(games.columns)
print(games.shape)

# Make a histogram of all the rating in the average_rating column.
plt.hist(games["average_rating"])
# Show the plot.
# plt.show()

# Print the first row of all the games with zero scores.
# The .iloc methods on the data frames allows us to index by position.
print(games[games["average_rating"] == 0].iloc[0])
# print the first row of all the games with scores greater than 0.
print(games[games["average_rating"] > 0].iloc[0])
# Remove any rows without user reviews.
games = games[games["users_rated"] > 0]
# Remove any rows with missing values.
games = games.dropna(axis=0)

# Clustering
# Initialized the model with 2 parameters -- number of clusters and random state (random seed).
kmeans_model = KMeans(n_clusters=5, random_state=1)
# Get only the numeric columns from games.
good_columns = games._get_numeric_data()
# Fit the model using the good columns.
kmeans_model.fit(good_columns)
# Get the cluster assignments.
labels = kmeans_model.labels_

# PCA, Plotting clustering.
# Create a PCA model.
pca_2 = PCA(2)
# Fit the PCA model on the numeric columns from earlier.
plot_columns = pca_2.fit_transform(good_columns)
# Make a scatter plot of each game, shaded according to cluster assignment.
plt.scatter(x=plot_columns[:,0], y=plot_columns[:,1], c= labels)
# show the plot.
# plt.show()

# Finding correlations
# Return pairwise correlation data frame
print(games.corr()["average_rating"])

# Picking predictor columns
# Get all the columns from the dataframe.
columns1 = games.columns.tolist()
# Filter the columns to remove ones we don't want.
columns1 = [c for c in columns1 if c not in ["bayes_average_rating", "average_rating", "type", "name"]]
# Store the variable we'll be predicting on.
target = "average_rating"
print(columns1)
print(target)

# Splitting into train and test sets
# Generate the training set.  Set random_state to be able to replicate results.
training = games.sample(frac=0.8, random_state=1)
# Select anything not in the training set and put it in the testing set.
testing = games.loc[~games.index.isin(training.index)]
# Print the shapes of both sets.
print(testing.shape)
# For train_test_split test
train, test = train_test_split(games, train_size=0.8, random_state=1)
print(test.shape, train.shape)

# Fitting a linear regression
# Initialize the model class.
model1 = LinearRegression()
model2 = LinearRegression()
# Fit the model to the training data
model1.fit(training[columns1], training[target])
model2.fit(train[columns1], train[target])

# Predicting error
# Generate out predictions for the test set.
predictions1 = model1.predict(testing[columns1])
predictions2 = model2.predict(test[columns1])
# Compute error between our test predictions and the actual values.
mse1 = mean_squared_error(predictions1, testing[target])
mse2 = mean_squared_error(predictions2, test[target])
print(mse1, mse2)

