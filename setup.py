from setuptools import setup

setup(
    name="enterprise_gwecc",
    version="0.1.0",
    description="Computes pulsar TOA delays due to gravitational waves from eccentric supermassive binary sources.",
    author="Abhimanyu Susobhanan",
    author_email="abhisrkckl@gmail.com",
    python_requires=">=3.8",
    package_dir={"": "src"},
    py_modules=["enterprise_gwecc"],
    install_requires=["numpy", "matplotlib", "enterprise-pulsar", "juliacall"],
)
