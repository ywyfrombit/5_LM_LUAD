# %%
import numpy as np 
import random 
from sklearn.datasets import load_iris  
from sklearn.model_selection import cross_val_score,GridSearchCV, KFold, train_test_split 
from sklearn.tree import DecisionTreeClassifier  
from deap import base, creator, tools, algorithms
import matplotlib.pyplot as plt
import warnings
from sklearn.exceptions import FitFailedWarning
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler,OneHotEncoder 
import pandas as pd
import matplotlib.pyplot as plt

#%matplotlib inline
from sksurv.datasets import load_breast_cancer,load_flchain, load_gbsg2
from sksurv.linear_model import CoxPHSurvivalAnalysis, CoxnetSurvivalAnalysis
from sklearn import set_config
from sklearn.impute import SimpleImputer
from sksurv.functions import StepFunction
from sksurv.metrics import (
    concordance_index_censored,
    concordance_index_ipcw,
    cumulative_dynamic_auc,
    integrated_brier_score,
)
from sksurv.nonparametric import kaplan_meier_estimator
from sksurv.util import Surv 
  
set_config(display="text")  # displays text representation of estimators

#%%

LUAD_data_file = "G:/MR_mGWAS/LUAD_RNASeq/LUAD_TPM_Clin.csv"
LUAD_data = pd.read_csv(LUAD_data_file ,header= 0)

LMn_file = "G:/LUAD_scRNASeq/Fig_and_chart/LM_gene_LUAD.csv"
LMn = pd.read_csv(LMn_file ,header= 0,index_col=0)
print(LMn)
LUAD_data = LUAD_data[~((LUAD_data["OS"]==0) & (LUAD_data["OS.time"]<30))]
LUAD_data = LUAD_data.dropna(subset=["OS.time"])
LM_gene =   LMn.index.tolist()
LUAD_cols = list(set(LUAD_data.columns).intersection(LM_gene))
X = LUAD_data.loc[:,LUAD_cols]
y_data = LUAD_data.loc[:,["OS","OS.time"]]
print(X.shape)
y_data["OS"] = y_data["OS"].apply(lambda x:True if x==1 else False)
y = y_data.to_records(index=False)
print(y)


#%%  

def evaluate_fitness(individual):  
    selected_features = [i for i, bit in enumerate(individual) if bit]  
    Xt = X.iloc[:, selected_features]  
     
    coxnet_pipe = make_pipeline(StandardScaler(), CoxnetSurvivalAnalysis(n_alphas=1, alphas=[0.2], alpha_min_ratio='auto',
                                                                     l1_ratio=0.9, max_iter=100))
    warnings.simplefilter("ignore", UserWarning)
    warnings.simplefilter("ignore", FitFailedWarning)
    coxnet_pipe.fit(Xt, y)
    estimated_alphas = coxnet_pipe.named_steps["coxnetsurvivalanalysis"].alphas_
    cv = KFold(n_splits=5, shuffle=True, random_state=0)
    gcv = GridSearchCV(
        make_pipeline(StandardScaler(), CoxnetSurvivalAnalysis(l1_ratio=0.5)),
        param_grid={"coxnetsurvivalanalysis__alphas": [[v] for v in estimated_alphas]},
        cv=cv,
        error_score=0.5,
        n_jobs=1,
    ).fit(Xt, y)

    cv_results = pd.DataFrame(gcv.cv_results_)
    score = cv_results.mean_test_score    
      
    return score  


# %%

creator.create("FitnessMax", base.Fitness, weights=(1.0,))  
creator.create("Individual", list, fitness=creator.FitnessMax)  
  
toolbox = base.Toolbox()  
toolbox.register("attr_bool", random.randint, 0, 1)  
toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_bool, n=X.shape[1])  
toolbox.register("population", tools.initRepeat, list, toolbox.individual)
toolbox.register("evaluate", evaluate_fitness)  
toolbox.register("mate", tools.cxTwoPoint)  
toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)  
toolbox.register("select", tools.selTournament, tournsize=3)  

stats = tools.Statistics(key = lambda ind: ind.fitness.values)
stats.register("avg",np.mean,axis=0)
stats.register("std",np.std,axis=0)
stats.register("min",np.min,axis=0)
stats.register("max",np.max,axis=0)
  

pop = toolbox.population(n=100)  
CXPB, MUTPB, NGEN = 0.5, 0.2, 200 
 
pop,log = algorithms.eaSimple(pop, toolbox, CXPB, MUTPB, NGEN, stats,verbose=False)  
  

best_ind = tools.selBest(pop, 1)[0]  
best_features = [i for i, bit in enumerate(best_ind) if bit]  
print(best_features)

log.header = "gen","nevals","avg","std","min","max"
print(log)

avgs = log.select("avg")
avgs_value = [arr.item() for arr in avgs]
print(avgs)


#%%
plt.plot(avgs_value,marker=None,linestyle="-")
plt.title("Fitness in GA iteration") 
plt.xlabel("Iteration") 
plt.ylabel("Fitness C-index")
plt.ylim(0.62,0.66) 
plt.grid(True)
plt.savefig('G:/LUAD_scRNASeq/GA-EnCox/fitness.pdf', format='pdf', dpi=300)
plt.show()


# %%

import pickle
with open("G:/LUAD_scRNASeq/GA-EnCox/result.pkl",'wb') as f:
    pickle.dump((best_features,avgs_value),f)


best_ind  = pd.DataFrame(best_ind)
best_ind.columns=["select_or_not"]
best_ind.to_csv("G:/LUAD_scRNASeq/GA-EnCox/best_features.csv",index=False)
# %%
