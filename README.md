
# Can Bitcoin Enhance Portfolio Performance?

This project examines whether Bitcoin and other cryptocurrencies improve portfolio diversification, inflation hedging, and risk-adjusted returns when combined with traditional financial assets.

## Overview

The analysis uses daily financial data from 2015 to 2025 for a diversified set of assets, including:

- Bitcoin, Ethereum, Dogecoin, and XRP
- Gold
- S&P 500
- 5-Year and 10-Year U.S. Treasury yield proxies
- Major foreign exchange pairs: EUR/USD, USD/JPY, GBP/USD, and AUD/USD

The project applies portfolio optimization techniques using rolling-window estimation and monthly rebalancing. It evaluates how cryptocurrencies behave in a multi-asset portfolio under inflationary and post-pandemic market conditions.

## Methods

The project is implemented in R and includes:

- Data collection from Yahoo Finance using `quantmod`
- CPI-based inflation adjustment using FRED data
- Daily return construction for price assets and Treasury yield proxies
- Rolling 252-day estimation window
- Monthly portfolio rebalancing
- Long-only portfolio constraints
- Trading-speed smoothing
- Portfolio performance evaluation

Four portfolio strategies are compared:

1. Equally Weighted Portfolio
2. Mean-Variance Portfolio
3. Maximum Sharpe Ratio Portfolio
4. Mean-CVaR Portfolio

## Key Features

- Inflation-adjusted return analysis
- Portfolio optimization using Modern Portfolio Theory
- Mean-CVaR downside risk optimization
- Sharpe Ratio and Adjusted Sharpe Ratio evaluation
- Certainty Equivalent return calculation
- Diversification metrics including Effective Number of Assets and Diversification Ratio
- Turnover and transaction cost estimation
- Cumulative wealth visualization

## Tools and Packages

The project uses the following R packages:

- `quantmod`
- `PerformanceAnalytics`
- `PortfolioAnalytics`
- `ROI`
- `ROI.plugin.quadprog`
- `ROI.plugin.glpk`
- `quadprog`
- `xts`
- `zoo`
- `tidyverse`
- `ggplot2`
- `flextable`
- `officer`

## Main Findings

The results show that cryptocurrencies can improve portfolio performance and diversification, but they do not consistently act as reliable inflation hedges. Gold remains a stronger traditional hedge, while Bitcoin and other digital assets behave more as return amplifiers within diversified portfolios.

The Equally Weighted portfolio performs strongly in terms of cumulative wealth, while optimization-based portfolios such as Mean-CVaR provide better downside risk control and lower drawdowns.

## Project Structure

```text
Capstone Final Code.R     # Main R script for data collection, optimization, and evaluation
README.md                 # Project documentation
````

## Author

Rabin Mishra<br/>
M.S. Economics<br/>
Texas A&M University

