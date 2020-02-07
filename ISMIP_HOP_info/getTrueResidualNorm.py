import re
import matplotlib.pyplot as plt
import numpy as np


def getTrueResidualNorm():  # Return all True residual norms printed in the text document 'textfile.txt'

    file = open("textfile.txt", "r")  # Text file with text output from running ISMIP_HOM_A_2d.py
    text = file.read()
    file.close()
    groups = re.findall("(true resid norm )(\d.\d+)(e\+\d{2})?", text)
    resid = []
    for i in groups:

        exponent = i[2]
        if len(exponent) < 3:  # Check if the exponent is on form e.g. 'e+05' or not
            resid.append(float(i[1]))
        else:
            resid.append(float(i[1]) * (10**int(exponent[3])))  # Multiplies with exponent
    return resid


def getLinearSolverAndPreconditioner():
    # Return type of linear solver and preconditioner printed in the text document 'textfile.txt'
    file = open("textfile.txt", "r")
    text = file.read()
    file.close()
    linear_solver = re.findall("(Linear solver used: )([a-z]+)", text)[0][1]
    preconditioner = re.findall("(Preconditioner used: )([a-z]+)", text)[0][1]
    return linear_solver, preconditioner


residuals = np.array(getTrueResidualNorm())
solver = getLinearSolverAndPreconditioner()[0]
preconditioner = getLinearSolverAndPreconditioner()[1]

x = np.arange(0,residuals.shape[0])

title_font = {'weight': 'bold', 'size': 20}  # Set fonts for the plots
x_label_font = {'size': 14}
y_label_font = {'size': 14}
x_tick_font = {'size': 12}
y_tick_font = {'size': 12}

figure = plt.figure(figsize=(12, 6), dpi=150)

plt.plot(x, residuals, label="solver: " + solver + ", preconditioner: " + preconditioner)
plt.title("True residual norm", title_font)
plt.xlabel("Krylov iteration", **x_label_font)
plt.ylabel("Residual norm", **y_label_font)
plt.legend()
plt.show()