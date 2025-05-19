import polars as pl
from blockhead.df_manip import read_vcf
from pathlib import Path

def test_read_csv_prints_correctly(capsys):
    path = Path(__file__).parent / "data" / "mier_test_1.vcf"
    df = read_vcf(str(path))
    assert isinstance(df, pl.DataFrame)