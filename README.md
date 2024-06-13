# MMED-PROJECT
The code repository for the MMED Project
# SEIAR Infectious Disease Model

## Introduction
This repository contains the code and data used for simulating the spread of an infectious disease using the SEIAR (Susceptible-Exposed-Infectious-Asymptomatic-Recovered) model. The model is applied to understand the impact of different symptomatic proportions on disease dynamics and vaccination strategies.

## Research Question
How do different proportions of symptomatic and asymptomatic individuals affect the spread of an infectious disease, and what are the implications for vaccination coverage needed to achieve herd immunity?

## System in Context
The SEIAR model is a compartmental model used in epidemiology to simulate the progression of individuals through different stages of infection and recovery. This study focuses on:
- Simulating the spread of disease for different proportions of symptomatic and asymptomatic individuals.
- Calculating the effective reproduction number ($R_{eff}$) and required vaccination coverage to control the disease.

## Methods
The model equations and parameters are implemented using R and the `deSolve` package. The differential equations are solved using the `lsoda` function. Simulations are conducted for various scenarios with different symptomatic proportions.

## Results
Key results include:
- Cumulative incidence of symptomatic infections for different scenarios.
- Required vaccination coverage to reduce $R_{eff}$ below 1 for each scenario.

## Discussion
The results highlight the critical role of symptomatic proportions in disease spread and vaccination strategies. Higher symptomatic proportions require higher vaccination coverage to achieve herd immunity.
