# OFP-Repair: Repairing Floating-point Errors via Original-Precision Arithmetic

This repository contains the data, code, and results for this work, as well as instructions on how to generate the results. This project is based on ACESO: https://bitbucket.org/FP-Aceso/aceso/src/master/

# Build the Docker Environment

Build the image from the repository root:

```bash
docker pull ubuntu:22.04
cd OFP-Repair
docker build -t ofp .
```

Run a container:

```bash
docker run -it --network=host --name ofp_container ofp /bin/bash
```

Inside the container, compile the project:

```bash
make
```

# Reproduce results

## RQ1

### Generate `condition_number_func_list.csv`

This file computes condition numbers with different step for all 63 functions.

```bash
python3 condition_number.py
```

## RQ2

### Generate `result.csv`

This file evaluates the original function, the Aceso patch, and the OFP-Repair
patches.

```bash
bin/testground.out fullinfo_ofp 10 10
```

The arguments are:

- `10`: Expansion maximum degree for OFP-Repair patches
- `10`: repeat count; metrics are averaged over 10 runs

You can also evaluate one function:

```bash
bin/testground.out fullinfo_ofp "1 exp_BI" 10 10 
```

### Generate `result_autornp.csv`

This file evaluates AutoRNP patches.

```bash
bin/testground.out fullinfo_autornp 10
```

The argument `10` is the repeat count. The output currently covers the AutoRNP
patches available in `autornp_patch.o`.

### Generate `result_herbie.csv`

This file evaluates Herbie patches from `herbie_optimized.c`.

```bash
bin/testground.out fullinfo_herbie 10
```

The argument `10` is the repeat count.

### Code for 14 cases unsupported by all baseline

See `/OFP-Repair/14cases_unsupported_by_all_baselines/`

## RQ3

### Generate `result_differencet_order.csv`

This file evaluates only OFP-Repair patches under different Expansion degrees.

```bash
bin/testground.out fullinfo_differencet_order 10 10
```

The arguments are:

- `10`: evaluate Expansion degrees from 1 through 10
- `10`: repeat count; metrics are averaged over 10 runs for each degree

## Case study

You the reproduce the fixing process for the three GSL bugs. 

```
python3 gsl_case_study.py
```

# Results

Please see `/OFP-Repair/Result/`, which contain the results of RQ1, RQ2, RQ3. 
