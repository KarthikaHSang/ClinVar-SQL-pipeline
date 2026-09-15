import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import classification_report, confusion_matrix

df = pd.read_csv("ml_features.csv")

# Encode categorical features as numbers, since scikit-learn needs numeric input
le_gene = LabelEncoder()
le_type = LabelEncoder()
le_review = LabelEncoder()

df["gene_encoded"] = le_gene.fit_transform(df["gene_symbol"])
df["type_encoded"] = le_type.fit_transform(df["type"])
df["review_encoded"] = le_review.fit_transform(df["review_status"])
df["in_domain_encoded"] = df["in_domain"].astype(int)

features = ["position", "gene_encoded", "type_encoded", "review_encoded",
            "number_submitters", "in_domain_encoded"]

X = df[features]
y = df["label"]

# Split into training and test sets - 80% to train on, 20% held back to
# evaluate on data the model has never seen
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42, stratify=y
)

print("Training set size:", len(X_train))
print("Test set size:", len(X_test))

model = RandomForestClassifier(n_estimators=100, random_state=42, class_weight="balanced")
model.fit(X_train, y_train)

y_pred = model.predict(X_test)

print("\nClassification report:")
print(classification_report(y_test, y_pred, target_names=["Benign", "Pathogenic"]))

print("Confusion matrix:")
print(confusion_matrix(y_test, y_pred))

print("\nFeature importances:")
importances = pd.Series(model.feature_importances_, index=features).sort_values(ascending=False)
print(importances)