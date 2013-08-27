====================
ADDING A NEW FEATURE
====================

1. Create an adapter class using script file in adapters/generators:

  cd adapters/generators
  bash runme.sh FeatureName "(int param1, int param2)"
  mv FeatureName.* ..

2. Add the new adapter to CMakeLists

3. Add the new adapter's header file to ConvertImageND.cxx

4. Add an invocation of the new adapter to ImageConverter::ProcessCommand

5. Edit the code in the adapters () operator

6. Update the documentation on the Wiki
