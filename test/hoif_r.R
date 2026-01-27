n <- 4000  # GPU 可以试 8000 / 10000
cat("Matrix size:", n, "x", n, "\n\n")

H1 <- rnorm(n)
H2 <- matrix(rnorm(n * n), n, n)
H3 <- rnorm(n)
tensors <- list(H1, H2, H2, H3)
expr <- "a,ab,bc,c->"

A <- sweep(H2, 1, H1, "*")
AA <- sweep(H2, 2, H3, "*")
diag(A) <- 0
diag(AA) <- 0
AAA <- A %*% AA
true <- (sum(AAA) - sum(diag(AAA)))/(n*(n-1)*(n-2))


# -------------------------------------------------
# Torch (auto → GPU if available)
# -------------------------------------------------
cat("Running Torch backend (auto dtype)...\n")
t1 <- system.time({
  res_torch_gpu <- ustat(
    tensors = tensors,
    expression = expr,
    backend = "torch",
    dtype = NULL   # 自动 float32(GPU) / float64(CPU)
  )
})
print(t1)
cat("Result (torch auto):", res_torch_gpu, "\n\n")
reticulate::py_last_error()  # 看完整 traceback，通常能看到 weight 来自哪里
# -------------------------------------------------
# Torch forced CPU float64
# -------------------------------------------------
cat("Running Torch backend (CPU float64)...\n")
t2 <- system.time({
  res_torch_cpu <- ustat(
    tensors = tensors,
    expression = expr,
    backend = "torch",
    dtype = "float64"
  )
})
print(t2)
cat("Result (torch cpu):", res_torch_cpu, "\n\n")

# -------------------------------------------------
# NumPy float64
# -------------------------------------------------
cat("Running NumPy backend (float64)...\n")
t3 <- system.time({
  res_numpy <- ustat(
    tensors = tensors,
    expression = expr,
    backend = "numpy",
    dtype = "float64"
  )
})
print(t3)
cat("Result (numpy):", res_numpy, "\n\n")

# -------------------------------------------------
# 数值差异对比
# -------------------------------------------------
cat("==== Result Differences ====\n")
cat("torch(GPU) vs true:", abs(res_torch_gpu - true), "\n")
cat("torch(CPU) vs true     :", abs(res_torch_cpu - true), "\n")
cat("torch(CPU) vs true     :", abs(res_torch_cpu - true), "\n")

cat("\n==== Speed Summary (seconds) ====\n")
print(rbind(
  torch_auto = t1[3],
  torch_cpu  = t2[3],
  numpy      = t3[3]
))

cat("\nDone.\n")
