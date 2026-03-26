# Dynamic Social Network Modeling: Relational Event Models & Hawkes Processes

This repository explores the statistical modeling of dynamic interactions within social networks. Rather than treating networks as static structures, this project focuses on the timing and sequence of events to understand the underlying drivers of human communication.

## 🎯 Project Overview
The core of this research is the application of **Relational Event Models (REM)**. The goal was to estimate the "rate" at which individuals interact based on their past history. 

I implemented a progression of three models, each adding a layer of social complexity:
1. **Model 1 (Accumulation):** A baseline model where the interaction rate depends solely on the total number of previous interactions between individuals.
2. **Model 2 (Recency):** An evolution of the first model that incorporates the influence of the most recent event, capturing short-term "reaction" patterns.
3. **Model 3 (Self-Excitement/Hawkes):** The most advanced stage, integrating a **Hawkes Process** into the REM framework. This captures the "bursty" nature of social interactions, where one event significantly increases the probability of another occurring shortly after.

## 💡 Key Methodologies
* **Point Processes:** Used to model interactions as events occurring in continuous time.
* **Self-Exciting Systems:** Applied Hawkes processes to quantify how the intensity of a relationship decays over time after an interaction.
* **Simulation & Validation:** Included scripts to simulate synthetic event sequences to validate the models' predictive accuracy.

## 🛠️ Tech Stack
* **Language:** R
* **Concepts:** Relational Event Models (REM), Hawkes Processes, Maximum Likelihood Estimation (MLE), Network Dynamics.

## 📂 Repository Structure
* **ASN_report.pdf**: The full technical report detailing the mathematical derivations, the three-model comparison, and the final results.
* **First model: only count.R**: Implementation of the frequency-based REM.
* **Second model: count with influence last event.R**: REM implementation with recency effects.
* **Third model: Hawkes in REM.R**: Integration of self-exciting point processes into the relational model.
* **simulate_hawkes.R**: A utility script for simulating event data following a Hawkes distribution.

## 📉 Results
The project demonstrates that social interactions are rarely uniform; they are driven by a combination of long-term history and short-term "bursts." Model 3, by incorporating the Hawkes decay, provided the most nuanced understanding of these temporal dynamics.
