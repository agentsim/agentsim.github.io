---
layout: post
title:  "The parallel evaluation of IIR filters - An Implementation"
date:   2020-09-16 13:49:17 -0400
---

## Introduction

This is an attempt to understand Raph Levien's [IIR filters can be evaluated in parallel](https://raphlinus.github.io/audio/2019/02/14/parallel-iir.html) by a person with no particular knowledge of DSP, filtering or relevant mathematics, just an interest in evaluating an IIR filter in parallel and a good deal of programming experience! The goal is to produce an implementation of this article demonstrating in code the approach Raph describes. The code won't be optimal by any means, but should serve as a sufficient base for building something much better.

Right from the get go we get an inkling this might not be as straight-forward as we might hope:

> Author’s note: I got slightly stuck writing this, so am publishing a somewhat unfinished draft

Throughout the article focuses on a simple RC-lowpass filter, given by the equation:

> `y[i] = x[i] * c + y[i - 1] * (1 - c)`

Raph also points out that this formulation of an IIR filter, while easy to understand has been superseded by the state space approach:

> Rather than writing it it in direct form, which emphasizes the serial evaluation strategy, it’s better to use matrices. Basically, y becomes a state vector, and a becomes a matrix. This is known as the state space approach to filters and has many advantages. I’d go so far to say that direct form is essentially obsolete now, optimizing for the number of multiplies at the expense of parallelism, numerical stability, and good modulation properties.

The state space approach will become important later. It is worth taking a little detour to clear it up.

---
## State Space Detour

All IIR filters are difference equations, which for our purposes are recursive equations that rely on one or more outputs from previous steps. Simply speaking, the deeper the recursion (more previous outputs), the more "powerful" it can be. This recursive level is the `order` of the filter. For example, to create a low-pass filter with a steeper cut-off, like a brick-wall filter, it will need more than just the single level of recursion of the RC-lowpass filter -- it will have a higher order.

A possible building block of such high-order IIRs is the biquad. A biquad is a second order (order = 2) filter. Biquads can be chained together to create even higher-order filters, a representation called second order sections.

This gets us to the `direct form` referenced in the quote above. Direct Form is a way to represent a biquad. One of the Direct Forms uses the following equation:

`y[n] = b0 * x[n] + b1 * x[n - 1] + b2 * x[n - 2] - a1 * y[n - 1] - a2 * y[n - 2]`

This is a second order filter because it relies on two historical values (`y[n - 1]` and `y[n - 2]`.) Thinking of the RC-lowpass filter in this way tells us that `b0` is `c` and `a1` is `c - 1` and the other constants are zero (our filter is first-order.)

### Conversion to State Space

All IIR filters can be representated by a transfer function, written as `H(z)`. [Wiki's](https://en.wikipedia.org/wiki/Digital_biquad_filter) introduction to biquads shows the transfer function. Mapping it to the direct form is trivial, so it is omitted here.

We want to convert the biquad form to a state space form. More generally it would be useful to be able to convert any IIR filter that can be written in terms of `H(z)` to state space form. This can be done by hand, but fortunately it can be done easily via GNU Octave's or Matlab's [tf2ss](https://www.mathworks.com/help/signal/ref/tf2ss.html) function:

```
c = 0.44
num = [c]
den = [1 c-1]
[A, B, C, D] = tf2ss(num, den)
A: 0.56
B: 1
C: 0.44
D: 0
```

Finally we can write the filter in state-space form. State-space form is generally written as:

```
x[n+1] = Ax[n] + Bu[n]
y[n] = Cx[n] + Du[n]
```

So far our equations have had `y` as the output (i.e the filtered signal) and `x` as the input (i.e. the signal to be filtered.) State-space changes this, the signal to be filtered is `u` and `x` is the _state vector_. Generally, `A` and `C` are matrices and `B` and `D` are vectors. Since our filter is first-order, each of these values are simply numbers. `D` is zero so it can be removed entirely and `B` is one so it can be simplified to:

```
x[n+1] = Ax[n] + u[n]
y[n] = Cx[n]
```

substituting the constant `C` from the RC-lowpass filter equation:

```
x[n+1] = (1 - C)x[n] + u[n]
y[n] = Cx[n]`
```

---

## Monoid Homomorphism Time

This is where the meat of the approach lies. Quoting from the article:

> The basic insight is that the target monoid is a function from input state to output state, and represents any integral number of samples. The monoid binary operator is function composition, which is by nature associative and has an identity (the identity function).
> In general, the amount of state required to represent such a function, as opposed to a single state value, is intractable, but in [this case] it works [...] because the function is linear.

This is going to require some unpacking.

---
### What is a monoid?

A monoid is a triplet of:
* a set being operated on
* an associative binary operator
* an identity element

The set might be the set of real numbers, for example. A binary operator is a function that takes two parameters that are elements of the set. Associativity is the property `(a op b) op c = a op (b op c)`, i.e. you can place the brackets arbitrarily. An identity element is an element `e` such that `e op x = x op e = x`. It is a no-op element.

For example, you might have a monoid representing addition: `(f64, +, 0)` meaning the set being operated on are the 64-bit floats, the binary operation is addition and the identity is `0` because `0+x = x+0 = x`. Similarly you could have a monoid representing multiplication: `(i32, *, 1)`.

---

Back to the article. Before tackling the homomorphism and what it means, we'll try to understand what Raph gives us as the path forward:

>The target monoid is a function of this form:
> `y_out = a * y_in + b`
>This can obviously be represented as two floats, we can write the representation as simply `(a, b)`. The homomorphism then binds a single sample of input into a function. Given the filter above, an input of `x` maps to `(c, x * (1 - c))`.
>Similarly, we can write out the effect of function composition in the two-floats representation space. Given `(a1, b1)` and `(a2, b2)`, their composition is `(a1 * a2, a2 * b1 + b2)`. Not especially complicated or difficult to compute.

We have to understand:

1. What is a monoid homomorphism?
2. How does the input mapping result in `(c, x * (1 - c)`?
3. Why is `y_out = a * y_in + b` our target monoid?
4. Why is our target monoid represented as `(a , b)`?
5. How is the binary operation (function composition) applied?

We'll do it backwards, from easiest to hardest to understand.

### How is the binary operation (function composition) applied?

Raph makes the following statement:

>Given `(a1, b1)` and `(a2, b2)`, their composition is `(a1 * a2, a2 * b1 + b2)`.

What is function composition? Say we have two functions, `f1` and `f2`, composition is simply `f2(f1(x))`. In our case the two functions are:
```
f1: y[n] = a1 * y[n-1] + b1
f2: y[n] = a2 * y[n-1] + b2
```

The function composition, `f2(f1(x))` is thus:
```
y[n] = a2 * (a1 * y[n-1] + b1) + b2
```
which multiplies out to:
```
y[n] = a2 * a1 * y[n-1] + a2 * b1 + b2
```
This is clearly still of the form `y[n] = a * y[n-1] + b`, with `(a, b)` now being `(a2 * a1, a2 * b1 + b2)`.


### Why is our target monoid represented as `(a , b)`?

This is now easier to understand. To represent the history of all prior invocations of the monoid, all we need are the `(a, b)` pair of the most recent function.

In other words, the history `fn(fn-1(fn-2(fn-3(...f0(0)...))))` can be represented as simply `fn(0)` provided we know the `(a, b)` pair.

Why do we know the initial input is zero? This chain of functions starts with `y[0] = ...`. But `y[0]` cannot be defined by a prior `y[-1]`; by definition there is nothing before the first element. Therefore the initial value for `y[n-1]` must be zero.

### Why is `y_out = a * y_in + b` our target monoid?

The state space representation for our filter was written as two equations:
```
x[n+1] = (1-C)x[n] + u[n]
y[n] = Cx[n]
```

Our goal is to evaluate the filter in parallel, so we want to get rid of recursion. The equation for the _state vector_, `x[n]` is recursive, so that's our target. We already know we can evaluate this equation without requiring the prior values provided we can compute values for `(a, b)`.

Also, from the introduction to this section, we're told:
>The basic insight is that the target monoid is a function from input state to output state

`x[n]` is the _state vector_ so our target monoid is evaluating the state of the filter at any `n` and producing the corresponding output state.

### How does the input mapping result in `(c, x * (1 - c)`?

All the pieces are in place, now we just need to compute those `(a, b)` values. We've got our target monoid:

`x[n+1] = (1 - C) * x[n] + u[n]`

let's write it as:

`y_out = (1 - C) * y_in + b`

Remember `u[n]` is the signal, so when we feed in a single sample `x`, it takes the place of `u[n]` which is `b` in our target form:

`y_out = (1 - C) * y_in + x`.

Pretty clearly, `(a, b)` is `(1 - C, x)`. That's not what we expected. There are two reasons for this.

* First, there's a mistake in the article. The mapping should result in `(1 - C, x * C)`. That's still not the result we have.
* Second, we are computing the state, not the final output. Once we have the state we still need to apply `y[n] = Cx[n]` to compute the final filtered result. Since that's not a recursive function it is trivial to compute in parallel. It turns out that `(1 - C, x * C)` will both compute the state and apply the `y[n]` equation. I'm not sure how to derive this result, I stumbled upon it by accident. Furthermore, in the general case I'm not sure `y[n]` can be computed in one step along with the state vector (consider the case where `D` is non-zero.)

### What is a Monoid Homomorphism?

A per usual, [Wiki](https://en.wikipedia.org/wiki/Monoid#Monoid_homomorphisms) is a good source for this.

In a nutshell, if you have two monoids with binary operators `M1` and `M2` and a function `f` that translates inputs for `M1` into `M2`, then you have a homomorphism if:

`f(M1(x, y)) = M2(f(x), f(y))` and given identity elements for the monoids `I1` and `I2` `f(I1) = I2`.

I'm not sure what the homomorphism is in this case since it requires two monoids and we've only discussed one.

---

## Prefix Sum Time

A prefix sum is a list of running totals, ex:
```
y[0] = x[0]
y[1] = x[0] + x[1]
y[2] = x[0] + x[1] + x[2]
```

This can be written as an equation that should seem pretty familiar by now:
`y[n] = y[n-1] + x[n]`

A prefix sum is not limited to addition, it applies to anything with a binary associative operator. That should also be familiar, any monoid can be the operator for a prefix sum.

It turns out that this can be computed in parallel. There are many algorithms for this. Here's an example of how it can be done (not at all efficient, but it illustrates that it is possible.) Let's say we have two threads and want to compute the prefix sum at each element for:

`[1, 2, 3, 4, 5, 6, 7, 8]`

We can divide this task into two:

The first thread computes: `[1, 3, 6, 10]`
The second computes: `[5, 11, 18, 26]`
Because addition is binary associative, we can now divide up the work of adding `10` to the second half:
The first thread computes: `[15, 21]`
The second thread computes: `[28, 36]`

This results in the completed prefix sum, computed in parallel with two threads:

`[1, 3, 6, 10, 15, 21, 28, 36]`

Now that we know a prefix sum can be computed in parallel, and that in general this can be done for any binary associative operator, and therefore any monoid, and since we know we can write our filter state equation as a monoid, then we know it can be computed in parallel.

---

## A Quick and Dirty Implementation

Let's start off with defining the RC-lowpass filter as a function. We'll evaluate this to test the output of the parallel implementation:

```
const C: f64 = 0.44;

fn rc_lowpass_filter(x: f64, y: f64) -> f64 {
    x * C + y * (1. - C)
}
```

Let feed this a simple 8-sample signal and see what we get:

```
fn main() {
    let input = [0.1, -0.2, 0.3, -0.4, 0.5, -0.6, 0.7, -0.8];
    let mut y_in = 0.;

    for (idx, x) in input.iter().enumerate() {
        y_in = rc_lowpass_filter(*x, y_in);
        println!("y[{}]: {}", idx, y_in);
    }
}
```
>```
>y[0]: 0.044000000000000004
>y[1]: -0.06336
>y[2]: 0.0965184
>y[3]: -0.12194969600000001
>y[4]: 0.15170817024
>y[5]: -0.17904342466560003
>y[6]: 0.20773568218726396
>y[7]: -0.2356680179751322
>```

Ok, now we have our sample input and sample output to test against.
Let's define our monoid state, a method to bind an input sample into a monoid state and our binary associative operator:

```
type State = (f64, f64);

fn bind(x: f64) -> State {
    (1. - C, x)
}

fn monoid_operator(f1: State, f2: State) -> State {
    let (a1, b1) = f1;
    let (a2, b2) = f2;

    (a1 * a2, a2 * b1 + b2)
}
```

All that should be pretty familiar by now. Now let's define a method to compute the prefix sum. This is taken from [Wiki](https://en.wikipedia.org/wiki/Prefix_sum#Algorithm_1:_Shorter_span,_more_parallel), it isn't the most efficient algorithm, but it serves our purpose. We'll use `into_par_iter()` from [Rayon](https://crates.io/crates/rayon) to get parallelism without much effort.

```
fn scan_monoid(input: &mut [f64]) {
    let log2n = (input.len() as f64).log2().floor() as usize;
    let mut x = vec![Vec::with_capacity(input.len()); log2n+1];

    for i in 0..input.len() {
        x[0].push(bind(input[i]));
    }

    for i in 0..log2n {
        let result = (0..input.len()).into_par_iter().map(|j| {
            if j < 2_usize.pow(i as u32) {
                x[i][j]
            } else {
                monoid_operator(x[i][j - 2_usize.pow(i as u32)], x[i][j])
            }
        }).collect();

        x[i + 1] = result;
    }

    for i in 0..input.len() {
        let (_a, b) = x[log2n][i];

        input[i] = C * b;
    }
}
```

Let's see what we get:

```
fn main() {
    let mut input = [0.1, -0.2, 0.3, -0.4, 0.5, -0.6, 0.7, -0.8];

    scan_monoid(&mut input);

    for (idx, y) in input.iter().enumerate() {
        println!("y[{}]: {}", idx, y);
    }
}
```
>```
>y[0]: 0.044000000000000004
>y[1]: -0.06336
>y[2]: 0.09651839999999999
>y[3]: -0.12194969600000001
>y[4]: 0.15170817024
>y[5]: -0.17904342466559997
>y[6]: 0.207735682187264
>y[7]: -0.23566801797513218
>```

Other than some floating point rounding differences, we get the same result.
