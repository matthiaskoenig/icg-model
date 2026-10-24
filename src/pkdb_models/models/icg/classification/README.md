# Classification Model for Survival after Hepatectomy

This classification quantifies the ability of the ICG-Model to predict survival in hepatectomy.
This is based on a prediction of postoperative ICG-R15 (and/or other pk-parameters in the future).

## Data in PK-DB: 
- icg/Seyama2009
- icg/Wakabayashi2004

## Levels of classification

1. - In clinical studies, a prediction of postoperative ICG-R15 is based on the liver remnant and the preoperative 
     ICG-R15: ICG-R15<sub>post-op</sub> = liver remnant * ICG-R15<sub>pre-op</sub>. A 2D classification model that takes
     into account the resection rate and the preoperative ICG-R15 is developed. 
     
     _Dataset_: (`survival` | `preop_r15` | `resection_rate`)
     
2. - A 1D classification model is developed, that only takes into account the postoperative ICG-R15 predicted as 
     described in(1).
     
     _Dataset_: (`survival` | `postop_r15_data`)
     
3. - With the ICG-Model, a prediction of postoperative ICG-R15 can be made. This prediction is calculated by mapping
     the preoperative ICG-R15 on the f_cirrhosis model-parameter (f_cirrhosis(ICG-R15<sub>pre-op</sub>)). With this 
     f_cirrhosis value and the planned resection rate postoperative ICG-R15 is predicted. A 2D classification model that
     takes into account the resection rate and the preoperative ICG-R15 (transformed into an f_cirrhosis-value) is
     developed. 
     
     _Dataset_: (`survival` | `f_cirrhosis` | `resection_rate`)
     
4. - A 1D classification model is developed, that only takes into account the postoperative ICG-R15 predicted by the 
     model as described in (2).
     
     _Dataset_: (`survival` | `postop_r15_model`)
    
These two classification models are compared. The key difference is the transformation of the preoperative ICG-R15 to 
an f_cirrhosis value in (2), which is where the model is used.

## Prepare data sets for classification

Create single dataset (TSV, pandas)
- `survival`: survival ([True|False])
- `study`: study identifier ([Seyama2009|Wakabayashi2004])
- `resection_rate`: resection rate ([0-1])
- `preop_r15`: preoperative ICG-R15 ([0-1])
- `f_cirrhosis`: f_cirrhosis ([0-1])  
- `postop_r15_data`: post-operative ICG-R15 ([0-1]) (calculated from resection rate & preop ICGR15 from data)
- `postop_r15_model`: post-operative ICG-R15 ([0-1]) (calculated from model)

## Classification Model
`scikit-learn` is used to train simple classifiers based on the clinical
data and ICG parameter predictions of the physiological-based pharmacokinetics model

see:
- https://scikit-learn.org/stable/auto_examples/classification/plot_classifier_comparison.html

**Evaluation metrics**
-> precision, recall, sensitivity, specificity, F1-score (balanced methods, e.g. balanced accuracy, f1-score instead of f-score)
-> matthews_corrcoef 
see:
- https://scikit-learn.org/stable/modules/model_evaluation.html#
- https://en.wikipedia.org/wiki/Confusion_matrixClassification
- https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-6413-7

Davis2006
    J. Davis, M. Goadrich, The Relationship Between Precision-Recall and ROC Curves, ICML 2006.
Flach2015
    P.A. Flach, M. Kull, Precision-Recall-Gain Curves: PR Analysis Done Right, NIPS 2015.


**Validation**
-> Bootstrap/K-Fold/Leave-one-out (Cross-validation):  
https://ogrisel.github.io/scikit-learn.org/sklearn-tutorial/modules/classes.html#module-sklearn.cross_validation
https://ogrisel.github.io/scikit-learn.org/sklearn-tutorial/modules/generated/sklearn.cross_validation.Bootstrap.html

**Trainings & performance curves**
https://scikit-learn.org/stable/modules/learning_curve.html

