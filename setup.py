from setuptools import setup

setup(
    name="enterprise_gwecc",
    version="0.1.1",
    description="Compute pulsar timing array signals due to eccentric supermassive binaries",
    author="Abhimanyu Susobhanan",
    author_email="abhisrkckl@gmail.com",
    python_requires=">=3.8",
    package_dir={"": "src"},
    py_modules=["enterprise_gwecc"],
    install_requires=["numpy", "scipy", "matplotlib", "enterprise-pulsar", "juliacall"],
)
