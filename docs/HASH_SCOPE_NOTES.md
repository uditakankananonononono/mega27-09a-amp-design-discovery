# Hash scope notes (2026-10-08 audit)

This sidecar note labels what the hash and checksum fields in the files below actually cover. The data files themselves are unchanged. It is a documentation clarification, not a license verdict and not a source re-admission.

Audited commit: `ec1b2c2c6cc04a6dd299566a443205db6de2f721`

## `data/raw/large_provenance.json`
- Lines (verified against the audited commit): 3, 7, 11, 15 (4 lines)
- Scope: The "sha256" values here are 16-hex-character prefixes, i.e. incomplete pins. They are not full SHA-256 verification.

