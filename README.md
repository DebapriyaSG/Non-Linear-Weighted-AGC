# Non-Linear-Weighted-AGC
`Enhancement_And_Quantitative_Evaluation_Proposed.m` enhances images using non-linear weight adjusted adaptive gamma correction.
`Enhancement_And_Quantitative_Evaluation_AGC.m` and `Quantitative_Evaluation_Original_Image.m` are added for comparison of enhancement.
In order to run, `Database` should be in the same directory where codes are downloaded.
`Enhancement_And_Quantitative_Evaluation_Medical_Images.m` enhances PET and MR slices (2D) using non-linear adaptive gamma correction.
NOTE: The `Medical_Data` directory is blank. This is because we do NOT yet have permission to share the data on which we worked.
In order to run `Enhancement_And_Quantitative_Evaluation_Medical_Images.m`, use MR and PET slices (2D) in dicom format.
The algorithm implemented here is published in `Multimedia Tools and Applications`, Volume 80, pp-3835-3862, 2021. In the published article, medical images
were not used.