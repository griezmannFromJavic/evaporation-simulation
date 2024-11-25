# Two-phase flow thermohydraulic oscillations
Continued work on [my master's thesis](https://zir.nsk.hr/islandora/object/fsb%3A9362).


# Governing Equations

## Mass Conservation Law
The mass conservation law is given by:

\[
\frac{\partial \rho}{\partial t} + \frac{\partial \rho u}{\partial x} = 0 \tag{1}
\]

---

## Momentum Conservation Law
The momentum conservation law is expressed as:

\[
\frac{\partial \rho u}{\partial t} + \frac{\partial \rho u^2}{\partial x} + \frac{\partial p}{\partial x} + \frac{\partial p}{\partial x}\bigg|_{\text{friction}} + \rho g \sin{\beta} = 0 \tag{2}
\]

---

## Energy Conservation Law
The energy conservation law can be written as:

\[
\frac{\partial}{\partial t}\left[\rho \left( h + \frac{u^2}{2} \right)\right] + \frac{\partial}{\partial x}\left[\rho u \left( h + \frac{u^2}{2} \right)\right] + \rho u g = \frac{\partial p}{\partial t} + q_w \tag{3}
\]

---

## Heat Conduction in the Pipe Wall
The equation describing the temperature of the pipe wall is the heat conduction equation, taken from \cite{galovic}. It describes heat conduction through a one-dimensional rod with a volumetric heat source term:

\[
\rho c \frac{\partial \theta_s}{\partial t} = - \frac{\partial q_x}{\partial x} + \phi_V
\]

\[
\frac{\partial \theta_s}{\partial t} = \frac{1}{\rho c} \frac{\partial}{\partial x} \lambda \frac{\partial \theta_s}{\partial x} + \frac{\phi_V}{\rho c}
\]

\[
\frac{\partial \theta_s}{\partial t} = \frac{\lambda}{\rho c} \frac{\partial^2 \theta_s}{\partial x^2} + \frac{\phi_V}{\rho c} \tag{4}
\]

