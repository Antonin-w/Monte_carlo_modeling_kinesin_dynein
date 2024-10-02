# Monte Carlo modeling of kinesis and dynein with Python

In this project, I've modeled kinesin and dynein movements on a microtubule using the Monte Carlo algorithm. The algorithm and values used have been taken from the following paper:

> Singh MP, Mallik R, Gross SP, Yu CC. [Monte Carlo modeling of single-molecule cytoplasmic dynein.](https://pubmed.ncbi.nlm.nih.gov/16103365/) Proc Natl Acad Sci U S A. 2005;102(34):12059-12064. doi:10.1073/pnas.0501570102

<br>

# Dynein 

![](res/dynein.gif)


```python
# Check dev/dynein.py to better understand functions created

# bind_unbind(s, ADP_released, F): function returning the new value of s depending if Binding/Unbinding occured. 

# hydrolysis_step(s, x): function returning new values of s, x and if ADP has been released or not, depending on Hydrolysis and if molecular motor took a step or no. 


for i in range (0, Nt-1):
    ADP_released = False
    F = Ktrap*x[i] # load
    t[i+1] = t[i] + delta_t
    s[i] = bind_unbind(s[i], ADP_released, F)

    if s[i] >= 1:
        s[i+1], x[i+1], ADP_released = hydrolysis_step(s[i], x[i])

    else:
        s[i+1], x[i+1], ADP_released = s[i], x[i], False
```

 <br /> 

# Kinesin 

![](res/kinesin.gif)

```python
for i in range (0, Nt-1):
    F = Ktrap*x[i] # load
    Pcat = (Kcat0*np.exp((-alpha*F*step)/(kB*T))) * delta_t # probability of catalysis of ATP

    t[i+1] = t[i] + delta_t

    # Binding/Unbinding ATP 
    p = uniform(0, 1)

    if s[i] == 0:
        s[i+1] = 1 if p <= Pon else 0 # Binding ATP with Pon probability
        x[i+1] = x[i] # MT doesn't move
    else:
        if p <= Pcat: # Hydrolisis takes place with probability Pcat
            epsilon = 1-(F/F0)**2
            p = uniform(0, 1)
            x[i+1] = x[i] + step if p <= epsilon else x[i] # Kinesis takes a step on the MT with probability epsilon
            s[i+1] = 0
        else: # Otherwise, stays at the same place
            s[i+1] = s[i]
            x[i+1] = x[i]

```

<br>

# Variables

```python
# Variables 

delta_t = 2e-4 # length of a time step
step = 8e-9 # taking 8-nm steps
Ktrap = 7e-6 # optical trap stiffness
F0 = 0.7e-11 # stalling force
Kcat0 = 55 # load-dependent rate constant
kB = 1.38064852e-23 # boltzmann constant
T = 300 # temperature (Kelvin)
alpha = 0.3 
beta = 0.7
ATP = 1e-3 # ATP concentration
d0 = 6e-9
Psyn0 = 0.23

# rate constants for binding
Kon1 = 4e5
Kon2 = 4e5
Kon3 = Kon2 / 4
Kon4 = Kon2 / 6

Koff1 = 10
Koff2 = 250
Koff3 = Koff2
Koff4 = Koff3

# probability of unbinding ATP
Poff1 = Koff1*delta_t
Poff2 = Koff2*delta_t
Poff3 = Koff3*delta_t
Poff4 = Koff4*delta_t
````


