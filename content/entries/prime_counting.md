---
title: A prime counting function with trigonometry
description: How could we construct an exact form for the prime counting
  function with simple  mathematical elements?
location: Oxford
date: 2024-01-01
draft: false
---
The **prime counting function (PCF)** counts the number of prime numbers less than or equal to a given number $x$. We denote it by $\Pi(x)$ and asymptotically obeys

$$
\Pi(x) \sim \frac{x}{\log(x)}, \quad x \to \infty.
$$

This is the **prime number theorem**, and it reveals the case that the density of primes decreases as we get to larger numbers. Computing this quantity exactly for arbitary $x$ is not straightforward. I here propose a ridiculous and intractable solution, which is exact in the limit

$$
\Pi(x; \epsilon, \rho) = \sum_{s=0}^{\rho} \frac{1}{1 + e^{2 \epsilon (s - x)}} \max \left[{\sigma(\epsilon,s), 0}\right],
$$
$$
\sigma(\epsilon,s) = \frac{\cos(\pi s)^{2 \epsilon}}{1 + e^{\epsilon (6 - 4 s)}} - \sum_{q=2}^{\rho} \frac{\cos\left(\frac{\pi s}{q}\right)^{2 \epsilon}}{1 + e^{\epsilon (6 q - 4 s)}}.
$$

In the positive infinite limit this actually tends to the true PCF

$$
\Pi(x) = \lim_{\rho, \epsilon \to \infty} \Pi(x; \epsilon, \rho), \quad \forall \rho, \epsilon \in \mathbb{N}.
$$

<!-- So how on earth does this have anything to do with primes? Let's consider the case of all numbers up until the prime number $11$. We will list these in the table below with all their divisors:

| 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |
|---|---|---|---|---|---|---|---|----|----|
| 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |
| 1 | 1 | 2 | 1 | 3 | 1 | 4 | 3 | 5  | 1  |
|   |   |   |   | 2 |   | 2 | 1 | 2  |    |
|   |   |   |   | 1 |   | 1 |   | 1  |    |

This is nice, but what it really tells us is that every number that is not a prime will eventually be divisible by another number that is smaller than it, excluding 1. Let's now place a point above a number line where its values are divisible by those we propose on a vertical axis. -->