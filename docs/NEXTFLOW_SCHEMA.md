# Understanding nextflow_schema.json

## What is it?

`nextflow_schema.json` is a **JSON Schema** file that defines all parameters for the scAnnex pipeline. It's part of the **nf-core standard** for Nextflow pipelines.

## Purpose

### 1. Parameter Documentation
Describes each parameter with:
- Name and type
- Default value
- Description
- Valid options (enums)
- Required vs optional

### 2. Input Validation
Nextflow automatically validates parameters against this schema:
```bash
# This will fail if min_genes is not a number
nextflow run main.nf --min_genes "invalid"
```

### 3. Automatic Help Generation
Powers the `--help` command:
```bash
nextflow run main.nf --help
```

Displays organized parameter groups with descriptions.

### 4. Integration with Tools

**Nextflow Tower / Seqera Platform:**
- Automatic GUI generation
- Parameter forms in web interface
- Preset configurations

**nf-core tools:**
- Schema linting
- Parameter documentation
- Website generation

## Structure

```json
{
  "definitions": {
    "input_output_options": {
      "properties": {
        "input": {
          "type": "string",
          "description": "Path to input file",
          "format": "file-path"
        }
      }
    },
    "quality_control_options": { ... },
    "normalization_options": { ... }
  }
}
```

## Should You Modify It?

**YES, when:**
- Adding new pipeline parameters
- Changing default values
- Adding parameter descriptions
- Updating valid options

**NO, don't:**
- Delete it (breaks nf-core compatibility)
- Remove structure (tools depend on it)

## Example: Adding a Parameter

If you add a new parameter to `nextflow.config`:

```groovy
params {
    new_feature = true
}
```

Update `nextflow_schema.json`:

```json
{
  "new_feature": {
    "type": "boolean",
    "default": true,
    "description": "Enable new feature",
    "fa_icon": "fas fa-star"
  }
}
```

## Location

**Must stay in project root** alongside `nextflow.config` and `main.nf`.

This is standard for all nf-core and nf-core-style pipelines.

## Tools

### Validate Schema
```bash
nf-core schema validate
```

### Build/Update Schema
```bash
nf-core schema build
```

### Lint Schema
```bash
nf-core schema lint
```

## For scAnnex

Your `nextflow_schema.json` is well-structured and defines:
- Input/output options
- Quality control parameters
- Normalization settings
- Integration options
- Clustering parameters
- Annotation settings

**Keep it in the root. It's important for:**
- Parameter validation
- Documentation generation
- Seqera Platform integration
- nf-core compatibility

## Summary

`nextflow_schema.json` is **essential infrastructure** for modern Nextflow pipelines. It provides structure, validation, and integration capabilities. Keep it in the root directory.
