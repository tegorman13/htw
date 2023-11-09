
import numpy as np
import pandas as pd
import random

def alm_response(input_value, c, input_layer, output_layer, weight_mat, train_vec=None):
    input_activation = np.exp(-c * (input_layer - input_value) ** 2)
    input_activation /= input_activation.sum()
    output_activation = weight_mat @ input_activation
    output_probability = output_activation / output_activation.sum()
    mean_response = (output_layer * output_probability).sum()
    return {
        'mean.response': mean_response, 
        'input.activation': input_activation, 
        'output.activation': output_activation
    }

def alm_update(cor_resp, c, lr, output_layer, input_activation, output_activation, weight_mat):
    fz = np.exp(-c * (output_layer - cor_resp) ** 2)
    teacher_signal = (fz - output_activation) * lr
    w_change = teacher_signal[:, np.newaxis] @ input_activation[np.newaxis, :]
    weight_mat += w_change
    weight_mat = np.clip(weight_mat, a_min=0, a_max=None)  # Replaces the commented out line
    return weight_mat

def alm_trial(input_value, cor_resp, c, lr, input_layer, output_layer, weight_mat):
    alm_resp = alm_response(input_value, c, input_layer, output_layer, weight_mat)
    updated_weight_mat = alm_update(
        cor_resp, c, lr, output_layer, 
        alm_resp['input.activation'], alm_resp['output.activation'], weight_mat
    )
    return {
        'mean.response': alm_resp['mean.response'], 
        'weight.mat': updated_weight_mat
    }

def alm_sim(dat, c, lr, input_layer, output_layer):
    weight_mat = np.full((len(output_layer), len(input_layer)), 0.000001)
    xt = dat['x']
    n = len(dat)
    st = np.zeros(n)  # Initialize the vector to store mean responses
    for i in range(n):
        trial = alm_trial(xt[i], dat['y'][i], c, lr, input_layer, output_layer, weight_mat)
        weight_mat = trial['weight.mat']
        st[i] = trial['mean.response']
    dat['almResp'] = st
    return {'d': dat, 'wm': weight_mat, 'c': c, 'lr': lr}

def exam_response(input_value, c, input_layer, output_layer, weight_mat, train_vec):
    nearest_train = train_vec[np.abs(input_value - train_vec).argmin()]
    aresp = alm_response(nearest_train, c, input_layer, output_layer, weight_mat)['mean.response']
    
    train_vec_sorted = np.sort(train_vec)
    nearest_index = np.where(train_vec_sorted == nearest_train)[0][0]
    x_under = train_vec_sorted[max(nearest_index - 1, 0)]
    x_over = train_vec_sorted[min(nearest_index + 1, len(train_vec_sorted) - 1)]
    
    m_under = alm_response(x_under, c, input_layer, output_layer, weight_mat)['mean.response']
    m_over = alm_response(x_over, c, input_layer, output_layer, weight_mat)['mean.response']
    
    if x_over == x_under:  # Prevent division by zero
        exam_output = aresp
    else:
        exam_output = round(aresp + ((m_over - m_under) / (x_over - x_under)) * (input_value - nearest_train), 3)
    
    return exam_output

def generate_data(x, type="linear", noise=np.nan):
    if type == "linear":
        y = np.round(2.2 * x + 30, 0)
    elif type == "exponential":
        y = np.round(200 * (1 - np.exp(-x / 25)), 0)
    elif type == "sinusoidal":
        y = np.sin(2 * np.pi * x)
    elif type == "quadratic":
        y = np.round(210 - ((x - 50) ** 2) / 12, 0)
    else:
        raise ValueError("type must be linear, exponential, quadratic, or sinusoidal")

    if not np.isnan(noise):
        y += np.round(np.random.normal(0, noise, len(y)), 2)

    return pd.DataFrame({'x': x, 'y': y, 'type': [type] * len(x)})

training_blocks = {
    'low': np.array([30.5, 36.0, 41.0, 46.5, 53.5, 59.0, 64.0, 69.5]),
    'med': np.linspace(30.0, 70.0, num=20),
    'high': np.linspace(30.0, 70.0, num=50)
}

def simulate_training_phase(function_type, density_level):
    # Generate training data based on the specified function and density level
    stimulus_magnitudes = training_blocks[density_level]
    training_data = generate_data(stimulus_magnitudes, type=function_type)
    
    # Repeat the stimuli for the required number of blocks
    block_repeats = 200 // len(stimulus_magnitudes)
    training_stimuli = np.tile(stimulus_magnitudes, block_repeats)
    
    # Shuffle the training stimuli for each block
    for i in range(block_repeats):
        start_index = i * len(stimulus_magnitudes)
        end_index = (i+1) * len(stimulus_magnitudes)
        np.random.shuffle(training_stimuli[start_index:end_index])

    return training_data

def simulate_transfer_phase():
    # Generate stimuli for the full range from 1 to 100 for the transfer phase
    all_stimuli = np.linspace(1.0, 100.0, num=100)
    transfer_data = pd.DataFrame({'x': all_stimuli})
    
    return transfer_data

def simulate_experiment(function_type, density_level):
    training_data = simulate_training_phase(function_type, density_level)
    transfer_data = simulate_transfer_phase()
    
    return training_data, transfer_data

# Example usage:
# training_data, transfer_data = simulate_experiment("linear", "low")
