import numpy as np
import numpy.linalg as LA
import matplotlib.pyplot as plt

def quadraticForm(_u):
    return \
        np.asscalar(
            np.transpose(_u) @ (A @ _u) + np.transpose(B) @ _u + C
        )


# Last sample of the long Markov Chain
u = np.matrix("-3.74175068e+00;-2.00791950e+01;-4.10018301e+01;9.84471817e+00;2.55054494e+01;2.95511739e+00;-4.18927875e+00;1.03589367e+00;"
              "1.11812336e+00;-1.44455201e+00;3.25228201e-01;-3.16556504e+00;9.69321155e-01;5.03833141e-01;-4.09041265e-01;-1.04627219e+06; "
              "1.04627201e+06")

u1 = u

# This is obtained in the Markov Chain as last point (belonging the sample above)
misfit = -1.51294100e+05

# Loading the matrices
A = np.loadtxt("A")
B = np.loadtxt("B")
C = np.loadtxt("C")

# If one has to check that the matrices have been loaded correctly.
# This is how they are used in the HMC sampler
print("A:\r\n", A)
print("B:\r\n", B)
print("C:\r\n", C, "\r\n\r\n")

print(np.linalg.eigvals(A))

# Calculate misfit to check
print("Last sample in long chain:")
print("misfit: ", quadraticForm(u))
print(u)
print("\r\n")

As = A + np.transpose(A)

rcond = 1e-32
u = -LA.pinv(As, rcond) @ B
u = np.expand_dims(u, 1)
u2 = u
print("Starting model:")
print("misfit: ", quadraticForm(u))
print(u)
print("\r\n")

# Split in positive/negative

# Plot parameters
plt.tight_layout()
x = np.linspace(1, 17, 17)
plt.scatter(x, (np.squeeze(np.asarray(u2))), marker='+', color="red")
plt.scatter(x, (np.squeeze(np.asarray(u1))), marker='+', color="blue")
plt.gca().set_yscale('symlog', linthreshy=0.01)
plt.legend(["Saddle point", "End of chain"])
plt.xlabel("Parameter number")
plt.ylabel("Parameter value (absolute)")
plt.gca().set_xticks(x)
plt.gca().grid(which='minor', alpha=0.2)
plt.gca().grid(which='major', alpha=0.5)
plt.show()
