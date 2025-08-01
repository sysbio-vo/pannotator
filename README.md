# A template for Nextflow pipeline projects

This template will provide a skeleton project for a Nextflow pipeline to be deployed to the farm, according to the
conventions developed by the PaM Informatics team.

[[_TOC_]]

## Template contents

### Nextflow Pipeline

#### Wrapper script

The `nextflow-template.sh` script provides a wrapper which will run the nextflow pipeline. This file must have the same
name as the git project (plus `.sh`).

#### "Hello world" nf pipeline.

The `main.nf` and `nextflow.config` files are included to provide a "hello world" pipeline.

#### Module template

A module template file `module.template` is provided. This includes placeholders of the form `{{placeholder_tag}}`
which are replaced in the deployed module file, according to the
[CI/CD pipeline configuration](https://gitlab.internal.sanger.ac.uk/sanger-pathogens/templates/ci-templates/-/blob/master/README.md?ref_type=heads#farm22pipeline-deployyml)

### CI/CD pipeline

The standard PaM Informatics nextflow deployment pipeline is included from the
[CI/CD templates repo](https://gitlab.internal.sanger.ac.uk/sanger-pathogens/templates/ci-templates/).

This includes checks to ensure that the pre-commit hooks have been run, some standard GitLab security scans, and
deployment on tagging.

#### Deployment

A test deployment will be done whenever a tag starting with `t` is created; a prod deployment is done whenever a
semantic version tag staring with `v` (e.g. v1.2.3) is created.

Note that the `.gitlab-ci.yml` file has a section that disables deployment, just to prevent the template itself from
being deployed to the farm. This section is marked by a comment, and should be deleted when the template is used for a
real project that is intended to be deployed.

Deployment includes:

- Copying the wrapper script to the farm, removing the `.sh`, and making it executable.
- Using the module template to create a module file on the farm.
- Copying all other files that do _not_ start with `.` to the farm; this will include the Nextflow files and any other scripts etc. that have been added to the project.

The deployment path is the value of the `FARM_PATH` variable in the deployment job, which should already be set to
the standard "custom install" path we use for dev/prod deployments.

### Workflow

Some standard files are provides to help with development workflow.

#### .gitignore file

This contains common patterns matching files that should be kept out of git. This can help to keep secrets out of the
repo, as well as clutter created by IDEs etc. Remember that secrets pushed to git by mistake will remain in the
history, even if deleted!

#### pre-commit configuration

A recommended pre-commit config is included (`.pre-commit.yml`).

Pre-commit hooks run before each commit automatically. We use them to auto-format code and to check for linting errors,
for example. W/o running the hooks, the CI pipeline may fail.

First, install pre-commit command itself if you don't have it already:

```
pip install pre-commit
```

Then, run the following from the repository's folder to install the pre-commit hooks:

```
pre-commit install
```

##### .talismanrc file

This file configures talisman, which is included in the pre-commit config. See the
[talisman documentation](https://github.com/thoughtworks/talisman/blob/main/README.md) to see how to deal with false
positives.

##### Running pre-commit hooks manually

[pre-commit provides various options](https://pre-commit.com/#pre-commit-run) for running hooks manually; for example,
so run all hooks on all files:

```
pre-commit run --all-files
```

## How to use this template

### Create a new project from the template

Start by creating a new project for your nextflow pipeline from the template

- In gitlab, click the blue "New project button"
- On the next screen, click on "Create from template"
- On the next screen, select the "Group" tab
- On the Group tab, next to "nextflow-template", click the blue "Use this template" button
- Fill in the form on the next screen

### Required changes in the project

The following manual changes _must_ be made in your new project

- Rename `nextflow-template.sh` to match your enw project name (plus `.sh` on the end).
- Edit the manifest in `nextflow.config` to replace the placeholder name, description etc.
- In `.gitlab-ci.yml`, delete the section which disables the deployment (see the comments in this file to identify the section you must delete).
- If you will be using python, uncomment the sections in `.pre-commit.yaml` for python-related hooks (see the comments in this file); and make corresponding changes in `.gitlab-ci.yml` to re-enable the jobs that check the python-related hooks (again, see the comments).
- If you want to make use of templated substitution strings, these can be specified (using sed scripts) in the [`./gitlab-ci/ci_file_sed_sub.tsv`](./gitlab-ci/ci_file_sed_sub.tsv) file, as documented for the `FARM_FILE_SED_SUB_FILE` [here](https://gitlab.internal.sanger.ac.uk/sanger-pathogens/templates/ci-templates/#options-for-all-deployments-pipeline_deploy_to_farm).

### Add your pipeline code

When the steps above have been completed, you should have a working project that will deploy a module to the farm (test
or prod, depending on the tag you use); but the pipeline provided will still be the "hello world" example.

You can now edit `main.nf` and `nextflow.config` to add a real pipeline, and add any additional scripts or other files
that are needed. The standard CI/CD piepline will copy all files that do not start with a `.` to the farm; it can
be configured easily to rename and change modes of file, and to make substitutions of text within files if required.
See the [CI/CD pipeline documentation](https://gitlab.internal.sanger.ac.uk/sanger-pathogens/templates/ci-templates/-/blob/master/README.md?ref_type=heads#farm22pipeline-deployyml)

Remember that if additional executables are provided, these should be added to the `module.template` file.
