import setuptools

# with open("README.md", "r") as fh:
#     long_description = fh.read()


def find_packages_with_dir(src_base, exclude):
    """Find packages under the given base directory and append their paths.
    """
    pkgs = setuptools.find_packages(src_base, exclude)
    return {pkg: src_base + '/' + pkg.replace('.', '/') for pkg in pkgs}


pkgs_dir = find_packages_with_dir('src/python', exclude=[])


setuptools.setup(
    name="chern",
    version="0.2",
    author="Ning Sun",
    author_email="ningsun.atom@gmail.com",
    description="Chern number calculator for 2D materials",
    # long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/atom-sun/chern.git",
    # packages=setuptools.find_packages(),
    packages=pkgs_dir.keys(),
    package_dir=pkgs_dir,
)
