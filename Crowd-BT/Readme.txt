=====================
 Crowd-BT Algorithm
=====================

1. General Information

This package implements the Crowd-BT algorithm for crowdsourced ranking inference.


2. Usage

To run a simulated experiment, open and run the file 'main.m'. You can change the first portion of 'main.m' (Basic inputs) for different experimental setups. The averaged accuracy will be saved.

Options of basic inputs include:

n_anno = number of workers
n_obj = number of items
budget = total available budget 
trials = number of independent trials to run

u, v = Beta(u,v) is the generating distribution of workers' reliability
alpha0, beta0 = Beta(alpha0,beta0) is the prior of workers' reliability

gamma = exploration-exploitaion tradeoff
