# Security Policy

## Supported versions

Only the latest release on the `main` branch receives fixes. Older
tagged releases are kept for reproducibility but are not patched.

## Reporting a vulnerability

If you believe you have found a security-relevant defect in PESTO,
please **do not open a public issue**. Instead, email the maintainer:

- Max Moldovan, `max.moldovan@adelaide.edu.au`

Include a minimal reproduction, the affected version, and the platform
(OS, R version, compiler). You will receive an acknowledgement within
five working days and, where applicable, a coordinated disclosure
timetable.

## Scope

In scope:

- Defects in PESTO’s R, C++, or shipped data that allow unintended file
  access, code execution, or numerical compromise of an inverse problem
  when the input is well-formed.
- Build-time defects in PESTO’s own sources that compromise installed
  artefacts.

Out of scope:

- Defects in upstream dependencies (Rcpp, Eigen, R itself). Report those
  to the relevant upstream project; we will coordinate where PESTO usage
  is the entry point.
- Behaviour that is documented as a constraint (e.g. examples that rely
  on user-provided file paths).
