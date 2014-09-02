#!/bin/bash

#make all the figures for paper on colors in group catalogues--Campbell 2014.

#figure 1
python mock_colors.py Mr19_age_distribution_matching_mock_sys_empty_shuffle_satrel_shuffle

#figure 2
#artwork made in keynote.

#figure 3
python mock_HOD_5.0.py Mr19_age_distribution_matching_mock_sys_empty_shuffle_satrel_shuffle

#figure 4
python mock_groups_f_L.py Mr19_age_distribution_matching_mock_sys_empty_shuffle_satrel_shuffle

#figure 5
python mock_groups_f_M.py Mr19_age_distribution_matching_mock_sys_empty_shuffle_satrel_shuffle

#figure 6
python mock_CLF_1.0.py Mr19_age_distribution_matching_mock_sys_empty_shuffle_satrel_shuffle

#figure 7
python conformity_w_errors_mock_runs_2.0.py Mr19_age_distribution_matching_mock Mr19_age_distribution_matching_mock_sys_empty_shuffle_satrel_shuffle

#figure 8
python HTP.py yang_mocks Mr19_age_distribution_matching_mock_sys_empty_shuffle_satrel_shuffle