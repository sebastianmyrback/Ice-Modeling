import numpy as np
import matplotlib.pyplot as plt


# Length of lists (number of iterations)

n1 = 12         # gmres, ilu
n2 = 11         # bicgstab, ilu
n3 = 10         # minres, none
n4 = 14         # tfqmr, ilu
n5 = 10         # tfqmr, jacobi
n6 = 9          # gmres, jacobi



number_of_solvers = 5



plt.figure(figsize=(12, 6), dpi=150)
title_font = {'weight': 'bold', 'size': 20}  # Set fonts for the plots
x_label_font = {'size': 14}
y_label_font = {'size': 14}

plt.title('Residual Norm', **title_font)
plt.xlabel('Number of Newton iterations', **x_label_font)
#plt.axis('tight')
plt.ylabel('Residual Norm', **y_label_font)
#file = open("newton_plot.txt", "r")
#print(file.read())



def plot_newtown():

    list_of_solvers = [0, '(gmres, ilu)', '(gmres, none)', '(bicgstab, ilu)', '(minres, none)', '(tfqmr, ilu)', '(tfqmr, jacobi)']
    num_of_iterations = [0, n1, n2, n3, n4, n5]
    list_of_iterations = []
    list_of_residuals = []
    counter = 0

    for n in range(len(num_of_iterations) - 1):

        if n == 0:
            current_iter = np.loadtxt("newton_plot.txt")[counter+1:counter + num_of_iterations[n + 1], 0]
            current_res = np.loadtxt("newton_plot.txt")[counter+1:counter + num_of_iterations[n + 1], 1]

        else:
            current_iter = np.loadtxt("newton_plot.txt")[counter+1:counter + num_of_iterations[n + 1], 0]
            current_res = np.loadtxt("newton_plot.txt")[counter+1:counter + num_of_iterations[n + 1], 1]

        print(current_iter)
        print(current_res)
        plt.plot(current_iter, current_res,
                 label=list_of_solvers[n + 1] + '; ' + str(num_of_iterations[n + 1]) + ' iterations')
        plt.legend(prop={'size': 12})
        counter += num_of_iterations[n + 1]
    plt.show()


def plot_special():
    list_of_solvers = [0, '(gmres, ilu)', '(bicgstab, ilu)', '(minres, none)', '(tfqmr, jacobi)']
    num_of_iterations = [0, n1, n2, n3, n5]
    list_of_iterations = []
    list_of_residuals = []
    counter = 0

    for n in range(len(num_of_iterations) - 1):

        if n == 0:
            current_iter = np.loadtxt("newton_plot_excl.txt")[counter+1:counter + num_of_iterations[n + 1], 0]
            current_res = np.loadtxt("newton_plot_excl.txt")[counter+1:counter + num_of_iterations[n + 1], 1]

        else:
            current_iter = np.loadtxt("newton_plot_excl.txt")[counter+1:counter + num_of_iterations[n + 1], 0]
            current_res = np.loadtxt("newton_plot_excl.txt")[counter+1:counter + num_of_iterations[n + 1], 1]

        print(current_iter)
        print(current_res)
        plt.plot(current_iter, current_res,
                 label=list_of_solvers[n + 1] + '; ' + str(num_of_iterations[n + 1]) + ' iterations')
        plt.legend(prop={'size': 12})
        counter += num_of_iterations[n + 1]
    plt.show()


plot_newtown()
#file.close()

