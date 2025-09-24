import polars as pl


class Clean:
    def __init__(self, file_path: str):
        # Read the Excel file and rename the first (unnamed) column to "protein"
        df = pl.read_excel(file_path)
        original_first = df.columns[0]
        self.df = df.rename({original_first: "protein"})
        print(f"Original shape: {self.df.shape}")

    def __step1__(self):
        """1) Drop rows with empty or null protein names."""
        before = self.df.height
        self.df = self.df.filter(
            pl.col("protein").is_not_null() & (pl.col("protein") != "")
        )
        after = self.df.height
        print(f"After step1 (empty names): {self.df.shape}, removed {before - after} rows")

    def __step2__(self):
        """2) Drop proteins whose all intensity values are 0 or null."""
        before = self.df.height
        value_cols = [c for c in self.df.columns if c != "protein"]
        mask = None
        for c in value_cols:
            cond = pl.col(c).is_null() | (pl.col(c) == 0)
            mask = cond if mask is None else (mask & cond)
        self.df = self.df.filter(~mask)
        after = self.df.height
        print(f"After step2 (all zeros/null): {self.df.shape}, removed {before - after} rows")

    def __step3__(self):
        """3) Drop contamination proteins (names containing 'contam')."""
        before = self.df.height
        self.df = self.df.filter(
            ~pl.col("protein").str.to_lowercase().str.contains("contam")
        )
        after = self.df.height
        print(f"After step3 (contaminants): {self.df.shape}, removed {before - after} rows")

    def __step4__(self):
        """4) Drop proteins with 'RepID=unknown'."""
        before = self.df.height
        self.df = self.df.filter(
            ~pl.col("protein").str.contains("RepID=unknown")
        )
        after = self.df.height
        print(f"After step4 (unknown RepID): {self.df.shape}, removed {before - after} rows")

    def __step5__(self):
        """5) Rename proteins starting with 'k77' to their RepID value."""
        before = self.df.height
        self.df = self.df.with_columns(
            pl.when(pl.col("protein").str.starts_with("k77"))
              .then(pl.col("protein").str.extract(r"RepID=([^|]+)", 1))
              .otherwise(pl.col("protein"))
              .alias("protein")
        )
        after = self.df.height
        print(f"After step5 (rename k77 RepIDs): {self.df.shape}, unchanged {before - after} rows")

    def __step6__(self):
        """6) Remove specific metagenome RepIDs."""
        blacklist = [
            "W1WC08_9ZZZZ", "W1WM01_9ZZZZ", "W1Y9K7_9ZZZZ", "W1YGV0_9ZZZZ",
            "W1YKP9_9ZZZZ", "W1YP67_9ZZZZ", "W1YRV8_9ZZZZ", "A0A0F9P276_9ZZZZ",
            "J9G8E8_9ZZZZ", "J9GD12_9ZZZZ"
        ]
        pattern = "|".join(blacklist)
        before = self.df.height
        self.df = self.df.filter(~pl.col("protein").str.contains(pattern))
        after = self.df.height
        print(f"After step6 (metagenome blacklist): {self.df.shape}, removed {before - after} rows")

    def __step7__(self):
        """7) Collapse duplicate proteins by summing intensities."""
        before = self.df.height
        value_cols = [c for c in self.df.columns if c != "protein"]
        self.df = (
            self.df
            .group_by("protein", maintain_order=True)
            .agg([pl.sum(c).alias(c) for c in value_cols])
        )
        after = self.df.height
        print(f"After step7 (collapse duplicates): {self.df.shape}, reduced from {before} to {after} rows")

    def run_filters(self, output_path: str = None):
        """Apply all cleaning steps in sequence."""
        for step in (
            self.__step1__, self.__step2__, self.__step3__, self.__step4__,
            self.__step5__, self.__step6__, self.__step7__
        ):
            step()
        print("✅ All filters applied. Final shape:", self.df.shape)
        if output_path:
            # Save cleaned dataframe to a CSV file
            # self.df = self.df.transpose()
            self.df.write_csv(output_path)
            print(f"Data saved to {output_path}")
        return self.df


if __name__ == "__main__":
    cleaner = Clean("test_all_patients_fragpipe_results_raw.xlsx")
    cleaned_df = cleaner.run_filters(output_path="CleanedTest.csv")
