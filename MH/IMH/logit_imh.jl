using DataFrames, GLM, Plots
temp = [53, 57, 58, 63, 66, 67, 67, 67, 68, 69, 70, 70, 70, 70, 72, 73, 75, 75, 76, 76, 78, 79, 81]
failure = [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0]

data = DataFrame(temp = temp, failure = failure)

logit_fit = glm(@formula(failure ~ temp), data, Binomial(), LogitLink())
plot(temp, predict(logit_fit), legend = false, xlabel = "Temperature", ylab = "Probability")
scatter!(temp, predict(logit_fit))