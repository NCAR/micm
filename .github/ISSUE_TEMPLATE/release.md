---
name: Version Release
about: Create an issue to make a new release
title: 'Release X.X.X'
labels: ''
assignees: ''

---

## Testing

- [ ] GitHub Actions are passing on `main`
- [ ] Update the version number in the include path in the README example
- [ ] Verify the README example compiles and runs successfully
- [ ] Confirm the README example output matches what is shown in the README

## Deployment

- [ ] Create a new branch (do **not** name it `release`)
- [ ] Update the version number in `CMakeLists.txt`
- [ ] Run the CMake configure step from the MICM source directory:
  ```
  cmake -B build -S .
  ```
  This forces CMake to write the correct version file into the include directory.
- [ ] Open `docs/sources/_static/switcher.json`
  - The first object points to the stable version — update its `name` to the new version number
- [ ] Update the version number in `CITATION.cff`
- [ ] On GitHub, merge `main` into `release` — **do NOT squash and merge**
  - Alternatively, merge locally and push: `git checkout release && git merge main && git push`
- [ ] Create a tag and add release notes on GitHub
  - Be sure to select the `release` branch as the target
