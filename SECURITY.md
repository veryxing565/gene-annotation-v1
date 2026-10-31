# Security and data handling

Use a trusted Linux environment and trusted local tools/images. Configuration,
FASTA/FASTQ, native source rules and external predictors are not untrusted-code
sandboxes. Do not run as root. Ownership/checksum checks reduce accidental
cleanup/corruption risk; they do not defend against a malicious local user.

Do not post private paths, reads, credentials or license keys in public issues.
For a security concern, contact the repository owner through an available
private GitHub channel; if none is available, open a minimal issue requesting
a private contact without publishing an exploit or private data.
