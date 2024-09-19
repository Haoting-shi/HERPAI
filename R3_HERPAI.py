import os
import argparse
import json
import pandas as pd
from autogluon.multimodal import MultiModalPredictor
from ray import tune
from autogluon.core.dataset import TabularDataset
from lifelines.utils import concordance_index

rst = {}

automm_hyperparameters = {
    "data.categorical.convert_to_text": False,
    "model.names": ["ft_transformer", "numerical_mlp"],
    "model.ft_transformer.embedding_arch": ["linear"],
    "env.batch_size": 128,
    "env.per_gpu_batch_size": 128,
    "env.eval_batch_size_ratio": 1,
    "env.num_workers": 12,
    "env.num_workers_evaluation": 12,
    "env.accelerator": "gpu",
    "env.num_gpus": 8,
    "optimization.max_epochs": 2000,  # Specify a large value to train until convergence
    "optimization.weight_decay": 1.0e-5,
    "optimization.lr_schedule": "polynomial_decay",
    "optimization.warmup_steps": 0.0,
    "optimization.patience": 20,
    "optimization.top_k": 3,
    "optimization.optim_type": 'adamw',
    "optimization.learning_rate": tune.uniform(0.00005, 0.005),
    "model.ft_transformer.num_blocks": tune.randint(3, 8),
    "model.ft_transformer.token_dim": tune.choice([96, 192, 256, 512]),
   # "model.ft_transformer.hidden_size": tune.choice([96, 192, 256, 512]),
   # "model.ft_transformer.ffn_hidden_size": tune.randint(4, 12),
    "model.ft_transformer.ffn_dropout": tune.uniform(0.0, 0.5),
    "model.ft_transformer.attention_n_heads": tune.choice([4, 8, 16, 32])

}

hyperparameter_tune_kwargs = {
    "searcher": "random",
    "scheduler": "FIFO",
    "num_trials": 80,
}

inputs = ['log_age',
          'log_bmi',
          'mense',
          'FH',
          'breast_surgery',
          'LN_surgery',
          'histology_group',
          'ER',
          'PR',
          'HER2_IHC',
          'Ki67',
          'T_stage',
          'N_stage',
          'grade_group']

for i in range(20):
    for j in range(5):
        impute_num = str(i+1)
        fold_num = str(j+1)

        train_path = f'./data/train/train_imputed{impute_num}_fold{fold_num}.csv'
        val_path = f'./data/validation/validation_imputed{impute_num}_fold{fold_num}.csv'
        test_path = f'./data/internal.test/internal.test_imputed{impute_num}.csv'

        train = TabularDataset(train_path)
        val = TabularDataset(val_path)
        test = TabularDataset(test_path)

        clean_train = train[inputs + ['iDFS']]
        clean_val = val[inputs + ['iDFS']]
        clean_test = test[inputs]

        predictor = MultiModalPredictor(
            label='iDFS',
            problem_type='regression',
            verbosity=1,
            path=f'./hyperparamtune_v5/model_impute{impute_num}_fold{fold_num}',
            enable_progress_bar=True,
            presets='best_quality_hpo'
        )

        ### model training
        predictor.fit(
            train_data=clean_train,
            tuning_data=clean_val,
            hyperparameters=automm_hyperparameters,
            hyperparameter_tune_kwargs=hyperparameter_tune_kwargs,
            save_path=f'./hyperparamtune_v5/fit_{impute_num}_{fold_num}',
            presets='best_quality',
        )

        pred = predictor.predict(data=clean_test, as_pandas=True)

        true_value = test['iDFS']

        c_index = concordance_index(true_value, pred)

        rst[str((impute_num, fold_num))] = c_index

        with open('./c-index-v5.json', 'w') as w:
            json.dump(rst, w)
