# Changelog

## 1.0.0 (2025-01-30)


### Features

* added authors, description, and reference section. closes [#2](https://github.com/MPUSP/snakemake-ms-proteomics/issues/2) ([c23716b](https://github.com/MPUSP/snakemake-ms-proteomics/commit/c23716bbe3de47906ae1d4a2bd2e0ebf3bfb09ed))
* added example test files ([722c721](https://github.com/MPUSP/snakemake-ms-proteomics/commit/722c721f1300ef8872adfdc7b4445f32d4ba58f5))
* added github actions ([83dcc44](https://github.com/MPUSP/snakemake-ms-proteomics/commit/83dcc44e94762f6fb4fe8a72852a2aee0c59557f))
* added log file export for all modules, closes [#1](https://github.com/MPUSP/snakemake-ms-proteomics/issues/1) ([94c1a5a](https://github.com/MPUSP/snakemake-ms-proteomics/commit/94c1a5a42264fdfe22f9a3b389043e99a9ce819a))
* added module for R MSstats package ([3a99ae1](https://github.com/MPUSP/snakemake-ms-proteomics/commit/3a99ae1047fee20b26aed0d180734d5a31dbecf7))
* added module to automatically retrieve proteome fasta file ([aa00ca9](https://github.com/MPUSP/snakemake-ms-proteomics/commit/aa00ca9b287ef7e9f50082bfcb6c690392fd9e77))
* added module to clean up temp files after pipeline execution; fixes [#1](https://github.com/MPUSP/snakemake-ms-proteomics/issues/1) ([785f877](https://github.com/MPUSP/snakemake-ms-proteomics/commit/785f877185bc4f92f650c8e8b4f87aa19cd901e8))
* added module to generate HTML report from R markdown template ([d239316](https://github.com/MPUSP/snakemake-ms-proteomics/commit/d239316a6a908036f0b4bb29bbc110fb128718f8))
* added module to sent email reports; closes [#8](https://github.com/MPUSP/snakemake-ms-proteomics/issues/8) ([04037a6](https://github.com/MPUSP/snakemake-ms-proteomics/commit/04037a6b71c0662efbf3588b188987f708dc5873))
* added more options to decoypyrat ([89db96c](https://github.com/MPUSP/snakemake-ms-proteomics/commit/89db96cdb53fc115b52e35dec1c2f425b405dc1d))
* added new module to check sample sheet ([b48cac4](https://github.com/MPUSP/snakemake-ms-proteomics/commit/b48cac43e99c8f31a61d811fddc47f745a0f1db8))
* added possibility to submit contrasts with samplesheet, fixes [#3](https://github.com/MPUSP/snakemake-ms-proteomics/issues/3) ([2f952e6](https://github.com/MPUSP/snakemake-ms-proteomics/commit/2f952e6bb3a0accdd41d70478aac454a204c768c))
* added protein summary tables to email report ([d9c2512](https://github.com/MPUSP/snakemake-ms-proteomics/commit/d9c2512baa62768a9b4aecc135f4214da1a73204))
* added rule to gather module and env log files; env info added to report ([f9b5897](https://github.com/MPUSP/snakemake-ms-proteomics/commit/f9b58979d0cf4237de5f5ff9635eee105adbd6a2))
* added schmematic overview of pipeline; closes [#2](https://github.com/MPUSP/snakemake-ms-proteomics/issues/2) ([df4aec0](https://github.com/MPUSP/snakemake-ms-proteomics/commit/df4aec0a5962b9631c05e72ab37beb43a1cd0066))
* added support for PDF export, closes [#7](https://github.com/MPUSP/snakemake-ms-proteomics/issues/7) ([ef5f75e](https://github.com/MPUSP/snakemake-ms-proteomics/commit/ef5f75efed6d9b274b38f8536da02c2586cd7dc1))
* added support to select workflow from sample type ([11dd081](https://github.com/MPUSP/snakemake-ms-proteomics/commit/11dd0817683a94f72f31bdeb10c8a3b7ded5d84b))
* added yml config for sm workflow catalog ([7cce96d](https://github.com/MPUSP/snakemake-ms-proteomics/commit/7cce96dd54fe4e01bee58d786fdb0497232533b1))
* change license to MIT ([d80702c](https://github.com/MPUSP/snakemake-ms-proteomics/commit/d80702cf01dee79acd8d145941a6b4fd71fc3dce))
* finished building basic HTML report with many different QC figures ([b9bd8c7](https://github.com/MPUSP/snakemake-ms-proteomics/commit/b9bd8c744cee4f0e3aba05ff896d213d2250c8ea))
* improved report, QC figures and documentation ([35c901a](https://github.com/MPUSP/snakemake-ms-proteomics/commit/35c901adf5aedd542aff39fd387ffedd16f3e0eb))
* major change in organisation, adding conda env definitions ([6def601](https://github.com/MPUSP/snakemake-ms-proteomics/commit/6def601ef53808eab0e1bb54d4ca878658b6fc1b))
* major update to repo structure and dependencies ([df35acc](https://github.com/MPUSP/snakemake-ms-proteomics/commit/df35acc8bf74c1773de6df397522098fc91d8cb4))
* streamlined snakefile regarding in/output paths, fixes [#6](https://github.com/MPUSP/snakemake-ms-proteomics/issues/6) ([c83e9f8](https://github.com/MPUSP/snakemake-ms-proteomics/commit/c83e9f8d995527703e67f4aa2fb26190c06d4e53))
* update catalog yml and snakemake options ([9198719](https://github.com/MPUSP/snakemake-ms-proteomics/commit/9198719f2bab44158ced00c82a3f7eb80f24ef5c))


### Bug Fixes

* corrected email address ([eae5ed2](https://github.com/MPUSP/snakemake-ms-proteomics/commit/eae5ed2bf9ce07fa1a5e95a4762bf0fe78b7132f))
* updated badges ([8094cfa](https://github.com/MPUSP/snakemake-ms-proteomics/commit/8094cfab7ef0285691f628f792b6d077c5f2a7f7))
