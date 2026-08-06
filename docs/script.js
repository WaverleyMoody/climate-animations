document.addEventListener('DOMContentLoaded', function () {
  const sections = document.querySelectorAll('.variable-section');

  sections.forEach(function (section) {
    const tabs = section.querySelectorAll('.tab');
    const panels = section.querySelectorAll('.panel');

    // Projection-specific links for this section only.
    const projectionLinks = section.querySelectorAll(
      '.code-link, .youtube-link'
    );

    tabs.forEach(function (tab) {
      tab.addEventListener('click', function () {
        const target = tab.getAttribute('data-projection');

        // Update this section's tabs only.
        tabs.forEach(function (t) {
          const isMatch = t === tab;

          t.classList.toggle('is-active', isMatch);
          t.setAttribute(
            'aria-selected',
            isMatch ? 'true' : 'false'
          );
        });

        // Show the selected projection panel.
        panels.forEach(function (panel) {
          const isMatch =
            panel.getAttribute('data-projection') === target;

          panel.classList.toggle('is-active', isMatch);

          // Pause videos in hidden panels.
          if (!isMatch) {
            const video = panel.querySelector('video');

            if (video && !video.paused) {
              video.pause();
            }
          }
        });

        // Update Python, R, and YouTube links.
        projectionLinks.forEach(function (link) {
          const url = link.getAttribute('data-' + target);

          if (url) {
            link.href = url;
          }
        });
      });
    });
  });

  // Playback speed controls.
  const speedButtons = document.querySelectorAll('.speed-btn');

  speedButtons.forEach(function (btn) {
    btn.addEventListener('click', function () {
      const panel = btn.closest('.panel');
      const video = panel
        ? panel.querySelector('video')
        : null;

      if (!video) return;

      video.playbackRate = parseFloat(
        btn.getAttribute('data-speed')
      );

      const siblingButtons =
        panel.querySelectorAll('.speed-btn');

      siblingButtons.forEach(function (b) {
        b.classList.toggle('is-active', b === btn);
      });
    });
  });

  // Loop controls.
  const loopButtons =
    document.querySelectorAll('[data-loop-toggle]');

  loopButtons.forEach(function (btn) {
    btn.addEventListener('click', function () {
      const panel = btn.closest('.panel');
      const video = panel
        ? panel.querySelector('video')
        : null;

      if (!video) return;

      video.loop = !video.loop;

      btn.classList.toggle(
        'is-active',
        video.loop
      );

      btn.setAttribute(
        'aria-pressed',
        video.loop ? 'true' : 'false'
      );
    });
  });

  // Read more / Show less controls.
  const aboutWraps =
    document.querySelectorAll('.meta-about-wrap');

  aboutWraps.forEach(function (wrap) {
    const text = wrap.querySelector('.meta-about');
    const toggle =
      wrap.querySelector('.meta-about-toggle');

    if (!text || !toggle) return;

    const isOverflowing =
      text.scrollHeight > text.clientHeight + 1;

    if (!isOverflowing) {
      toggle.style.display = 'none';
      return;
    }

    toggle.addEventListener('click', function () {
      const expanded =
        text.classList.toggle('is-expanded');

      toggle.textContent = expanded
        ? 'Show less'
        : 'Read more';
    });
  });
});
