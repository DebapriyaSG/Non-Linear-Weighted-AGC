# Non-Linear-Weighted-AGC
`Enhancement_And_Quantitative_Evaluation_Proposed.m` enhances images using non-linear weight adjusted adaptive gamma correction.
`Enhancement_And_Quantitative_Evaluation_AGC.m` and `Quantitative_Evaluation_Original_Image.m` are added for comparison of enhancement.
In order to run, `Database` should be in the same directory where codes are downloaded. \
`Enhancement_And_Quantitative_Evaluation_Medical_Images.m` enhances PET and MR slices (2D) using non-linear adaptive gamma correction.\
NOTE: The `Medical_Data` directory, mentioned in `Enhancement_And_Quantitative_Evaluation_Medical_Images.m` does not contain the required data. This is because we do NOT yet have permission to share the data on which we worked.
In order to run `Enhancement_And_Quantitative_Evaluation_Medical_Images.m`, use MR and/or PET slices (2D) in dicom format and store it in `Medical_Data` directory.
The algorithm implemented here is published in `Non-linear weight adjustment in adaptive gamma correction for image contrast enhancement. `_Multimedia Tools and Applications,_` Volume 80, pp 3835-3862, 2021`. In the published article, medical images
were not used.
