x = collect(range(0.0; length=301, stop=2.0 * pi))
y = sin.(x) + randn(length(x))
ns1 = NSplineBasis(x; df=5, intercept=false)

X = ns1(x)
fit1 = lm(X, y)

d = (; x=x, y=y)
fit2 = lm(@formula(y ~ 0 + ns(x, 5)), d)

newx = collect(0.0:0.5:3.5)
@test isapprox(predict(fit1, ns1(newx)), predict(fit2, (; x=newx)))
