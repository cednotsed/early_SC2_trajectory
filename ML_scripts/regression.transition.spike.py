# Import Statements
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import shap
from sklearn.model_selection import GridSearchCV, RepeatedKFold, cross_validate, KFold, cross_val_score
from sklearn.metrics import make_scorer, explained_variance_score, r2_score, mean_absolute_error
from xgboost import XGBRegressor

# Get paths
cwd = Path.cwd()
print(cwd)
datasets = cwd / '../results/prediction_out/mutation_stats/'
results = cwd / '../results/ML_out'

## Data Preprocessing
### Load data
dataset = sys.argv[1]
df = pd.read_csv(datasets / f'{dataset}.stats.csv')

print(dataset)

X = df.loc[:, ['n', 'median_freq', 'max_freq',
               'blosum62_score', 'delta_charge', 'abs_charge',
               'delta_mw', 'abs_mw', 'delta_hydropathy',
               'abs_hydropathy', 'mean_escape', 'delta_bind',
               'delta_expr']]

X['n'] = np.log10(X['n'])
X.mean_escape = X.mean_escape.fillna(-1)
X.delta_bind = X.delta_bind.fillna(-100)
X.delta_expr = X.delta_bind.fillna(-100)

y = np.log10(df.loc[:, 'global_n'] + 1)

## Model training and evaluation
def optimise_evaluate(X, y):
    np.random.seed(66)

    # Hyperparemeter Optimisation using grid search (F1)
    regressor = XGBRegressor()
    n_estimators = range(100, 1000, 200)
    max_depth = range(1, 10, 1)
    colsample_bytree = np.linspace(0.1, 1, 10)

    param_grid = dict(max_depth=max_depth,
                      n_estimators=n_estimators,
                      colsample_bytree=colsample_bytree,
                      n_jobs=[1])

    inner_cv = KFold(n_splits=10, shuffle=True)
    outer_cv = KFold(n_splits=10, shuffle=True)

    # Inner CV
    model = GridSearchCV(regressor,
                         param_grid,
                         scoring='neg_mean_squared_error',
                         n_jobs=6,
                         cv=inner_cv,
                         verbose=1)

    model.fit(X, y)
    best_params = model.best_params_
    print(best_params)

    # Custom metrics
    expl_var = make_scorer(explained_variance_score)
    rsquared = make_scorer(r2_score)
    mae = make_scorer(mean_absolute_error)

    scoring = {'r2': rsquared,
               'mae': mae}

    # Outer CV
    outer_results = cross_validate(model, X=X, y=y, cv=outer_cv, scoring=scoring, n_jobs=4)
    outer_results = pd.DataFrame(outer_results)

    return outer_results, best_params


# Tune hyperparameters
raw_results, raw_params = optimise_evaluate(X, y)

# Results
res = pd.DataFrame(raw_results).mean()[['test_r2', 'test_mae']]

# Fit final model
np.random.seed(66)
raw_model = XGBRegressor(**raw_params, enable_categorical=True)
raw_model.fit(X, y)

# Get prediction errors
y_pred = raw_model.predict(X)

test_results = pd.DataFrame({'y_test': y, 'y_pred': y_pred, 'mutation_name': df.mutation_name})
test_results.to_csv(results / f'within_dataset_results/{dataset}.spike.results.csv', index=False)

# SHAP analysis
X1000 = shap.utils.sample(X, 1000)

explainer_ebm = shap.Explainer(raw_model.predict, X1000)
shap_values_ebm = explainer_ebm(X)

# Parse SHAP values
data_df = pd.DataFrame(shap_values_ebm.data, columns=shap_values_ebm.feature_names)
shap_df = pd.DataFrame(shap_values_ebm.values, columns=shap_values_ebm.feature_names).add_suffix('.shap')
merged = pd.concat([df.loc[:, 'mutation_name'], data_df, shap_df], axis=1)

# Beeswarm plot
shap.plots.beeswarm(shap_values_ebm, show=False)

# Export results
shap.plots.beeswarm(shap_values_ebm, show=False)
merged.to_csv(results / f'shap_out/{dataset}.spike.shap.csv', index=0)
raw_results.to_csv(results / f'results_out/{dataset}.spike.results.csv', index=0)
plt.savefig(results / f'beeswarm_out/{dataset}.spike.pdf', bbox_inches='tight')
