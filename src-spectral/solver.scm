(define (stable-polynomial-evaluation a x)
  "Evaluate the polynomial u_N(x) = a_0 T_0(x) + ... + a_N T_N(x)
   using a stable recursive algorithm with Chebyshev polynomials."
  (let loop ((k (- (length a) 1))
             (B_N 0)
             (B_N_plus_1 0))
    (if (< k 0)
        ;; Final calculation: u_N(x) = a_0 + x * B_N - B_N_plus_1
        (+ (car a) (* x B_N) (- B_N_plus_1))
        ;; Recursive calculation: B_k = a_k + 2 * x * B_N - B_N_plus_1
        (loop (- k 1)
              (+ (list-ref a k) (* 2 x B_N) (- B_N_plus_1))
              B_N))))

